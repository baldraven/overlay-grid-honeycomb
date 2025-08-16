pub(crate) mod model;
pub(crate) mod routines;

use std::collections::VecDeque;

use crate::{
    model::{GeoVertices, Geometry2, IsIrregular, RefinementLevel, SiblingDartId},
    routines::{compute_overlapping_grid_size, detect_orientation_issue},
};
use honeycomb_core::{
    cmap::{CMap2, CMapBuilder, NULL_DART_ID},
    geometry::{CoordsFloat, Vertex2},
};
use thiserror::Error;
use vtkio::Vtk;

#[derive(Error, Debug)]
/// Errors that can occur during overlay grid operations
pub enum OverlayGridError {
    /// An orientation issue has been detected in the input geometry.
    #[error("boundary isn't consistently oriented - {0}")]
    InconsistentOrientation(&'static str),
    /// The specified geometry does not match one (or more) requirements of the algorithm.
    #[error("input shape isn't conform to requirements - {0}")]
    InvalidShape(&'static str),
    /// The VTK file used to try to build a `Geometry2` object contains invalid data
    /// (per VTK's specification).
    #[error("invalid/corrupted data in the vtk file - {0}")]
    BadVtkData(&'static str),
    /// The VTK file used to try to build a `Geometry2` object contains valid but unsupported data.
    #[error("unsupported data in the vtk file - {0}")]
    UnsupportedVtkData(&'static str),
}

fn dart_origin<T: CoordsFloat>(map: &CMap2<T>, dart: u32) -> Vertex2<T> {
    assert_ne!(dart, NULL_DART_ID);

    map.force_read_vertex(map.vertex_id(dart)).unwrap()
}

/// Collects all darts in a face by following the beta1 relation
fn collect_face_darts<T: CoordsFloat>(map: &CMap2<T>, start_dart: u32) -> Vec<u32> {
    let mut face_darts = Vec::new();

    // We want sure that our structure starts with a regular dart
    let mut current_dart = start_dart;
    while !is_regular(map, current_dart) {
        current_dart = map.beta::<1>(current_dart);
    }

    // then after we retrieve only regular darts
    let new_start = current_dart;
    loop {
        face_darts.push(current_dart);
        current_dart = canonical_beta1(map, current_dart).1;
        if current_dart == new_start {
            break;
        }
    }
    face_darts
}

fn canonical_beta1<T: CoordsFloat>(map: &CMap2<T>, dart: u32) -> (u32, u32) {
    assert_ne!(dart, NULL_DART_ID);

    let mut steps = 1;
    let mut current_dart = map.beta::<1>(dart);

    while !is_regular(map, current_dart) {
        steps += 1;
        current_dart = map.beta::<1>(current_dart);
    }

    (steps, current_dart)
}

fn is_regular<T: CoordsFloat>(map: &CMap2<T>, dart: u32) -> bool {
    assert_ne!(dart, NULL_DART_ID);

    let attr = map.force_read_attribute::<IsIrregular>(dart);
    attr.is_none()
}

/// Splits a boundary edge by inserting a new vertex at the midpoint
/// If the edge is lightly irregular, it doesn't do anything
fn split_boundary_edge<T: CoordsFloat>(
    map: &mut CMap2<T>,
    edge_dart: u32,
) -> Result<u32, OverlayGridError> {
    //let next = map.beta::<1>(edge_dart);
    let (n_irregular, next) = canonical_beta1(map, edge_dart);
    let dart1;

    if n_irregular == 2 {
        dart1 = map.beta::<1>(edge_dart);
        map.force_remove_attribute::<IsIrregular>(dart1);
        map.force_unlink::<1>(edge_dart).unwrap();
    } else if n_irregular == 1 {
        let opposite = map.beta::<2>(edge_dart);
        dart1 = if opposite != NULL_DART_ID {
            let dart1 = map.allocate_used_darts(2);
            let dart2 = dart1 + 1;
            map.force_write_attribute::<IsIrregular>(dart2, IsIrregular(true));
            let opposite_next = map.beta::<1>(opposite);

            map.force_unlink::<2>(edge_dart).unwrap();
            map.force_unlink::<1>(opposite).unwrap();
            map.force_link::<2>(dart1, opposite).unwrap();
            map.force_link::<2>(edge_dart, dart2).unwrap();
            map.force_link::<1>(opposite, dart2).unwrap();
            map.force_link::<1>(dart2, opposite_next).unwrap();
            dart1
        } else {
            map.allocate_used_darts(1)
        };

        map.force_unlink::<1>(edge_dart).unwrap();
        map.force_link::<1>(dart1, next).unwrap();

        let subdivde_point =
            Vertex2::<T>::average(&dart_origin(map, edge_dart), &dart_origin(map, next));
        map.force_write_vertex(dart1, subdivde_point);
    } else {
        //assert!(n_irregular == 4, "E: too much irregular");
        let sub_dart1 = map.beta::<1>(edge_dart);
        dart1 = map.beta::<1>(sub_dart1);
        map.force_remove_attribute::<IsIrregular>(dart1);
        map.force_unlink::<1>(sub_dart1).unwrap();
    }

    assert_ne!(dart1, NULL_DART_ID, "E: dart1 is NULL_DART_ID");
    Ok(dart1)
}

/// Creates inner edges from boundary to center
fn create_inner_edge<T: CoordsFloat>(
    map: &mut CMap2<T>,
    edge_dart: u32,
) -> Result<(u32, u32), OverlayGridError> {
    let going_to_center = map.allocate_used_darts(2);
    let going_from_center = going_to_center + 1;

    let edge_dart_next = map.beta::<1>(edge_dart);
    if edge_dart_next != NULL_DART_ID && !is_regular(map, edge_dart_next) {
        map.force_link::<1>(edge_dart_next, going_to_center)
            .unwrap();
    } else {
        map.force_link::<1>(edge_dart, going_to_center).unwrap();
    }
    map.force_link::<1>(going_to_center, going_from_center)
        .unwrap();

    Ok((going_to_center, going_from_center))
}

/// State management for subdivision process
#[derive(Default)]
struct SubdivisionState {
    going_to_center_prev: u32,
    going_from_center_first: u32,
    dart1_prev: u32,
}

impl SubdivisionState {
    fn update_for_iteration(&mut self, going_to_center: u32, dart1: u32) {
        self.going_to_center_prev = going_to_center;
        self.dart1_prev = dart1;
    }
}

/// Distributes geo vertices into 4 quadrants by sorting them in place
/// Returns the (start_index, length) tuple for each quadrant
fn distribute_geo_vertices_to_quadrants<T: CoordsFloat>(
    geo_verts: &mut [Vertex2<T>],
    parent_range: (u32, u32),
    midpoint: &Vertex2<T>,
) -> [(u32, u32); 4] {
    let start = parent_range.0 as usize;
    let end = start + parent_range.1 as usize;

    if start >= geo_verts.len() || end > geo_verts.len() {
        return [(0, 0), (0, 0), (0, 0), (0, 0)];
    }

    // Extract the slice we're working with
    let slice = &mut geo_verts[start..end];

    // Sort vertices into quadrants: bottom-left, bottom-right, top-right, top-left
    slice.sort_by(|a, b| {
        let a_quadrant = get_quadrant(a, midpoint);
        let b_quadrant = get_quadrant(b, midpoint);
        a_quadrant.cmp(&b_quadrant)
    });

    // Count vertices in each quadrant
    let mut quadrant_counts = [0u32; 4];
    for vertex in slice.iter() {
        let quadrant = get_quadrant(vertex, midpoint);
        quadrant_counts[quadrant] += 1;
    }

    // Calculate ranges for each quadrant
    let mut current_start = start as u32;
    let mut quadrant_ranges = [(0u32, 0u32); 4];

    for i in 0..4 {
        quadrant_ranges[i] = (current_start, quadrant_counts[i]);
        current_start += quadrant_counts[i];
    }

    quadrant_ranges
}

/// Determines which quadrant a vertex belongs to relative to a midpoint
/// Returns: 0=bottom-left, 1=bottom-right, 2=top-right, 3=top-left
fn get_quadrant<T: CoordsFloat>(vertex: &Vertex2<T>, midpoint: &Vertex2<T>) -> usize {
    let x_right = vertex.x() >= midpoint.x();
    let y_top = vertex.y() >= midpoint.y();

    match (x_right, y_top) {
        (false, false) => 0, // bottom-left
        (true, false) => 1,  // bottom-right
        (true, true) => 2,   // top-right
        (false, true) => 3,  // top-left
    }
}

// Fetches the 4 siblings dart to call refine_cell on the 4 of them
fn refine_with_pairing<T: CoordsFloat>(
    map: &mut CMap2<T>,
    working_dart: u32,
    geo_verts: &mut [Vertex2<T>],
    balance_pile: &mut VecDeque<u32>,
) -> Result<Vec<Vec<u32>>, OverlayGridError> {
    let face1 = map.face_id(working_dart);
    let face2 = map.force_read_attribute::<SiblingDartId>(face1).unwrap().0;
    let face3 = map.force_read_attribute::<SiblingDartId>(face2).unwrap().0;
    let face4 = map.force_read_attribute::<SiblingDartId>(face3).unwrap().0;

    let children1 = refine_cell(map, face1, geo_verts, balance_pile).unwrap();
    let children2 = refine_cell(map, face2, geo_verts, balance_pile).unwrap();
    let children3 = refine_cell(map, face3, geo_verts, balance_pile).unwrap();
    let children4 = refine_cell(map, face4, geo_verts, balance_pile).unwrap();

    Ok(vec![children1, children2, children3, children4])
}

/// Refines a single cell and returns the child cells that need further refinement
/// to understand better what needs to be done, see README.md
fn refine_cell<T: CoordsFloat>(
    map: &mut CMap2<T>,
    working_dart: u32,
    geo_verts: &mut [Vertex2<T>],
    balance_pile: &mut VecDeque<u32>,
) -> Result<Vec<u32>, OverlayGridError> {
    let face_id = map.face_id(working_dart);

    let face_darts = collect_face_darts(map, working_dart);
    assert_eq!(face_darts.len(), 4, "E: face darts count is not 4");

    // Check if one of the opposites already has a lower refinement level
    // If so we add it to the balance pile
    let current_depth = map
        .force_read_attribute::<RefinementLevel>(face_id)
        .unwrap()
        .0;
    for &dart in &face_darts {
        let opposite = map.beta::<2>(dart);
        if opposite == NULL_DART_ID {
            continue;
        }
        let opposite_depth = map
            .force_read_attribute::<RefinementLevel>(map.face_id(opposite))
            .unwrap()
            .0;

        if current_depth > opposite_depth {
            balance_pile.push_back(opposite);
        }
    }

    // Get the current cell's geo vertices range
    let geo_range = map.force_read_attribute::<GeoVertices>(face_id);
    let parent_range = if let Some(range) = geo_range {
        range.0
    } else {
        // No geo vertices associated with this face
        (0, 0)
    };

    let mut state = SubdivisionState::default();

    // Calculate cell bounds for quadrant sorting
    let cell_min = dart_origin(map, face_darts[0]);
    let cell_max = dart_origin(map, face_darts[2]);
    let midpoint = Vertex2::<T>::average(&cell_min, &cell_max);

    // Extract and sort geo vertices into quadrants if we have any
    let quadrant_ranges = if parent_range.1 > 0 {
        distribute_geo_vertices_to_quadrants(geo_verts, parent_range, &midpoint)
    } else {
        // No vertices to distribute
        [(0, 0), (0, 0), (0, 0), (0, 0)]
    };

    // Perform the cell subdivision
    for (i, &edge_dart) in face_darts.iter().enumerate() {
        // Boundary edge splitting
        let dart1 = split_boundary_edge(map, edge_dart).unwrap();

        // Inner edge splitting
        let (going_to_center, going_from_center) = create_inner_edge(map, edge_dart).unwrap();

        if i > 0 {
            // middle iteration
            map.force_link::<2>(going_from_center, state.going_to_center_prev)
                .unwrap();
            map.force_link::<1>(going_from_center, state.dart1_prev)
                .unwrap();
        }
        if i == 3 {
            // last iteration
            map.force_link::<2>(going_to_center, state.going_from_center_first)
                .unwrap();
            map.force_link::<1>(state.going_from_center_first, dart1)
                .unwrap();
        } else if i == 0 {
            // first iteration
            state.going_from_center_first = going_from_center;
            map.force_write_vertex(going_from_center, midpoint);
        }
        state.update_for_iteration(going_to_center, dart1);
    }

    // After subdivision, assign geo vertices to the 4 new faces and collect those that need refinement

    let mut child_cells_to_refine = Vec::new();

    // Assign geo vertices to each quadrant and track which cells need refinement
    for (i, &dart) in face_darts.iter().enumerate() {
        let child_face_id = map.face_id(dart);
        if quadrant_ranges[i].1 > 0 {
            // This child cell has geo vertices, so it will need refinement
            map.force_write_attribute::<GeoVertices>(
                child_face_id,
                GeoVertices(quadrant_ranges[i]),
            );
            child_cells_to_refine.push(dart);
        }
    }

    // Add the sibling attribute and the refinement level
    let refinement_level = map
        .force_read_attribute::<RefinementLevel>(face_id)
        .unwrap()
        .0;
    for (i, &dart) in face_darts.iter().enumerate() {
        let child_face_id = map.face_id(dart);
        map.force_write_attribute::<SiblingDartId>(
            child_face_id,
            SiblingDartId(map.face_id(face_darts[(i + 1) % 4])),
        );
        map.force_write_attribute::<RefinementLevel>(
            child_face_id,
            RefinementLevel(refinement_level + 1u16),
        );
    }

    // Return child cells that need further refinement
    Ok(child_cells_to_refine)
}

#[allow(clippy::needless_pass_by_value)]
/// Create an overlay grid from a VTK file input
pub fn overlay_grid<T: CoordsFloat>(
    file_path: impl AsRef<std::path::Path>,
    _grid_cell_sizes: [T; 2],
) -> Result<CMap2<T>, OverlayGridError> {
    // --- IMPORT VTK INPUT
    let geometry_vtk = match Vtk::import(file_path) {
        Ok(vtk) => vtk,
        Err(e) => panic!("E: could not open specified vtk file - {e}"),
    };
    //----/

    // --- BUILD OUR MODEL FROM THE VTK IMPORT
    let geometry = Geometry2::<T>::try_from(geometry_vtk)?;
    //----/

    // --- FIRST DETECTION OF ORIENTATION ISSUES
    detect_orientation_issue(&geometry)?;
    //----/

    // --- Overlapping grid size computation
    let [max_x, max_y, min_x, min_y] = compute_overlapping_grid_size(&geometry)?;
    let _quadtree_bounds = [min_x, min_y, max_x, max_y];
    const _MAX_DEPTH: usize = 5;
    let width = max_x - min_x;
    let height = max_y - min_y;

    let padding_factor = T::from(0.1).unwrap();
    let x_padding = width * padding_factor;
    let y_padding = height * padding_factor;

    let grid_min_x = min_x - x_padding;
    let grid_max_x = max_x + x_padding;
    let grid_min_y = min_y - y_padding;
    let grid_max_y = max_y + y_padding;

    /*     let mut geo_verts = vec![
           Vertex2::<T>(T::from(0.56).unwrap(), T::from(0.48).unwrap()),
           Vertex2::<T>(T::from(3.68).unwrap(), T::from(0.40).unwrap()),
           Vertex2::<T>(T::from(0.40).unwrap(), T::from(3.52).unwrap()),
           Vertex2::<T>(T::from(3.60).unwrap(), T::from(3.60).unwrap()),
           Vertex2::<T>(T::from(2.08).unwrap(), T::from(1.92).unwrap()),
           Vertex2::<T>(T::from(1.12).unwrap(), T::from(2.88).unwrap()),
           Vertex2::<T>(T::from(3.12).unwrap(), T::from(0.88).unwrap()),
           //Vertex2::<T>(T::from(0.01).unwrap(), T::from(0.01).unwrap()),
       ];
    */
    let mut geo_verts = geometry.vertices;

    // --- BUILD THE CMap2

    let mut map: CMap2<T> = CMapBuilder::<2, _>::from_n_darts(4)
        .add_attribute::<IsIrregular>()
        .add_attribute::<GeoVertices>()
        .add_attribute::<SiblingDartId>()
        .add_attribute::<RefinementLevel>()
        .build()
        .unwrap();

    map.force_link::<1>(1, 2).unwrap();
    map.force_link::<1>(2, 3).unwrap();
    map.force_link::<1>(3, 4).unwrap();
    map.force_link::<1>(4, 1).unwrap();

    // Set vertices based on geometry bounds with padding
    map.force_write_vertex(1, Vertex2(grid_min_x, grid_min_y)); // Bottom-left
    map.force_write_vertex(2, Vertex2(grid_max_x, grid_min_y)); // Bottom-right
    map.force_write_vertex(3, Vertex2(grid_max_x, grid_max_y)); // Top-right
    map.force_write_vertex(4, Vertex2(grid_min_x, grid_max_y)); // Top-left

    /*     map.force_write_vertex(1, Vertex2(T::from(0.0).unwrap(), T::from(0.0).unwrap()));
    map.force_write_vertex(2, Vertex2(T::from(4.0).unwrap(), T::from(0.0).unwrap()));
    map.force_write_vertex(3, Vertex2(T::from(4.0).unwrap(), T::from(4.0).unwrap()));
    map.force_write_vertex(4, Vertex2(T::from(0.0).unwrap(), T::from(4.0).unwrap())); */

    // putting all the geo vertices in the initial quadrant.
    map.force_write_attribute::<GeoVertices>(1, GeoVertices((0u32, geo_verts.len() as u32)));
    map.force_write_attribute::<RefinementLevel>(1, RefinementLevel(0));

    // Perform multi-level refinement
    const MAX_DEPTH: u32 = 3;
    geo_refinement(&mut map, &mut geo_verts, MAX_DEPTH)?;

    Ok(map)
}

/// Performs refinement of cells that contain geo vertices
/// Processes one level at a time to ensure proper ordering
fn geo_refinement<T: CoordsFloat>(
    map: &mut CMap2<T>,
    geo_verts: &mut [Vertex2<T>],
    max_depth: u32,
) -> Result<(), OverlayGridError> {
    let mut balance_pile = VecDeque::new();

    let to_refine = refine_cell(map, 1, geo_verts, &mut balance_pile).unwrap();

    let mut current_level = VecDeque::new();
    for dart in to_refine {
        current_level.push_back(dart);
    }

    println!("Starting refinement with 1 initial cell");

    // Process each refinement level
    for depth in 0..max_depth {
        if current_level.is_empty() {
            println!("No more cells to refine at depth {}", depth);
            break;
        }

        println!(
            "Processing refinement level {} with {} cells",
            depth,
            current_level.len()
        );

        let mut next_level = VecDeque::new();

        // Refine all cells at current level
        while let Some(working_dart) = current_level.pop_front() {
            let parent_cells =
                refine_with_pairing(map, working_dart, geo_verts, &mut balance_pile).unwrap();

            // Add children to next level
            for child_cells in parent_cells {
                for child_dart in child_cells {
                    next_level.push_back(child_dart);
                }
            }
        }

        println!(
            "Level {} complete. {} cells added for next level",
            depth,
            next_level.len()
        );

        balancing(map, &mut balance_pile, geo_verts);

        current_level = next_level;
    }

    if !current_level.is_empty() {
        println!(
            "Refinement stopped at max depth {} with {} cells remaining",
            max_depth,
            current_level.len()
        );
    }

    Ok(())
}

fn balancing<T: CoordsFloat>(
    map: &mut CMap2<T>,
    balance_pile: &mut VecDeque<u32>,
    geo_verts: &mut [Vertex2<T>],
) {
    while !balance_pile.is_empty() {
        let dart = balance_pile.pop_front().unwrap();
        refine_with_pairing(map, dart, geo_verts, balance_pile).unwrap();
    }
}

#[cfg(test)]
mod test;
