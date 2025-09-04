use std::collections::VecDeque;

use crate::model::{GeoVertices, IsIrregular, RefinementLevel, SiblingDartId};
use honeycomb_core::{
    cmap::{CMap2, NULL_DART_ID},
    geometry::{CoordsFloat, Vertex2},
};

use crate::routines::helpers::*;

/// Splits a boundary edge by inserting a new vertex at the midpoint
/// see README.md for visualisation
fn split_boundary_edge<T: CoordsFloat>(map: &mut CMap2<T>, edge_dart: u32) -> u32 {
    let (n_irregular, next) = canonical_beta1(map, edge_dart);
    // If we have an edge that has already been refined because of a refinement of its opposite
    // we need to re-use the extra darts and connect accordingly
    let dart1 = match n_irregular {
        4 => split_boundary_edge_n_4(map, edge_dart),
        2 => split_boundary_edge_n_2(map, edge_dart),
        1 => split_boundary_edge_n_1(map, edge_dart, next),
        _ => panic!("E: weirdly irregular ; a cell has probably been over-refined"),
    };
    assert_ne!(dart1, NULL_DART_ID, "E: dart1 is NULL_DART_ID");

    dart1
}

fn split_boundary_edge_n_1<T: CoordsFloat>(map: &mut CMap2<T>, edge_dart: u32, next: u32) -> u32 {
    // If we have an opposite we have to also split it and make its new connection, this is how we handle T junctions
    let opposite = map.beta::<2>(edge_dart);
    let dart1 = if opposite != NULL_DART_ID {
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

    dart1
}

fn split_boundary_edge_n_4<T: CoordsFloat>(map: &mut CMap2<T>, edge_dart: u32) -> u32 {
    let sub_dart1 = map.beta::<1>(edge_dart);
    let dart1 = map.beta::<1>(sub_dart1);

    // we'll use the pre-existing dart1, so it's not irregular anymore
    map.force_remove_attribute::<IsIrregular>(dart1);
    map.force_unlink::<1>(sub_dart1).unwrap();

    dart1
}

fn split_boundary_edge_n_2<T: CoordsFloat>(map: &mut CMap2<T>, edge_dart: u32) -> u32 {
    let dart1 = map.beta::<1>(edge_dart);

    // we'll use the pre-existing dart1, so it's not irregular anymore
    map.force_remove_attribute::<IsIrregular>(dart1);
    map.force_unlink::<1>(edge_dart).unwrap();

    dart1
}

/// Creates inner darts from boundary to center
fn create_inner_darts<T: CoordsFloat>(map: &mut CMap2<T>, edge_dart: u32) -> (u32, u32) {
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

    (going_to_center, going_from_center)
}

/// Fetches the 4 siblings dart to call refine_cell on the 4 of them
pub fn refine_with_pairing<T: CoordsFloat>(
    map: &mut CMap2<T>,
    working_dart: u32,
    geo_verts: &mut [Vertex2<T>],
    balance_pile: &mut VecDeque<u32>,
) -> Vec<u32> {
    let face1 = map.face_id(working_dart);
    let face2 = map.force_read_attribute::<SiblingDartId>(face1).unwrap().0;
    let face3 = map.force_read_attribute::<SiblingDartId>(face2).unwrap().0;
    let face4 = map.force_read_attribute::<SiblingDartId>(face3).unwrap().0;

    let children1 = refine_cell(map, face1, geo_verts, balance_pile);
    let children2 = refine_cell(map, face2, geo_verts, balance_pile);
    let children3 = refine_cell(map, face3, geo_verts, balance_pile);
    let children4 = refine_cell(map, face4, geo_verts, balance_pile);

    vec![children1, children2, children3, children4]
}

/// Refines a single cell and returns the child cells that need further refinement
/// to understand better what needs to be done, see README.md
pub fn refine_cell<T: CoordsFloat>(
    map: &mut CMap2<T>,
    working_dart: u32,
    geo_verts: &mut [Vertex2<T>],
    balance_pile: &mut VecDeque<u32>,
) -> u32 {
    let face_id = map.face_id(working_dart);
    let face_darts = collect_face_darts(map, working_dart);
    assert_eq!(face_darts.len(), 4, "E: face darts count is not 4");

    // Step 1: Update balance pile with neighbors that need balancing
    update_balance_pile_for_neighbors(map, face_id, &face_darts, balance_pile);

    // Step 2: Prepare geo vertices distribution
    let parent_geo_range = get_parent_geo_range(map, face_id);
    let quadrant_ranges =
        prepare_geo_vertices_for_subdivision(map, geo_verts, &face_darts, parent_geo_range);

    // Step 3: Perform the actual cell subdivision
    perform_cell_subdivision(map, &face_darts);

    // Step 4: Update child cell attributes and find cells needing further refinement

    update_child_cell_attributes(map, &face_darts, &quadrant_ranges, face_id)
}

/// Updates the balance pile with neighboring cells that need balancing
fn update_balance_pile_for_neighbors<T: CoordsFloat>(
    map: &CMap2<T>,
    face_id: u32,
    face_darts: &[u32],
    balance_pile: &mut VecDeque<u32>,
) {
    let current_depth = map
        .force_read_attribute::<RefinementLevel>(face_id)
        .unwrap()
        .0;

    for &dart in face_darts {
        let opposite = map.beta::<2>(dart);
        if opposite == NULL_DART_ID {
            continue;
        }

        let opposite_depth = map
            .force_read_attribute::<RefinementLevel>(map.face_id(opposite))
            .unwrap()
            .0;

        if current_depth > opposite_depth {
            add_to_balance_pile_if_not_exists(map, opposite, balance_pile);
        }
    }
}

/// Adds a face to the balance pile if no sibling is already present
fn add_to_balance_pile_if_not_exists<T: CoordsFloat>(
    map: &CMap2<T>,
    opposite_dart: u32,
    balance_pile: &mut VecDeque<u32>,
) {
    // Check if any sibling of this opposite face is already in the balance pile
    // TODO: find a more efficient way of doing this
    // because we don't want to add 2 siblings in the balance pile
    let opposite_face_id = map.face_id(opposite_dart);

    if !is_any_sibling_in_balance_pile(map, opposite_face_id, balance_pile) {
        balance_pile.push_back(opposite_dart);
    }
}

/// Checks if any of the 4 siblings is already in the balance pile
fn is_any_sibling_in_balance_pile<T: CoordsFloat>(
    map: &CMap2<T>,
    face_id: u32,
    balance_pile: &VecDeque<u32>,
) -> bool {
    let mut current_sibling = face_id;

    for _ in 0..4 {
        if balance_pile
            .iter()
            .any(|&dart| map.face_id(dart) == current_sibling)
        {
            return true;
        }
        current_sibling = map
            .force_read_attribute::<SiblingDartId>(current_sibling)
            .unwrap()
            .0;
    }

    false
}

/// Gets the geo vertices range for the current cell
fn get_parent_geo_range<T: CoordsFloat>(map: &CMap2<T>, face_id: u32) -> (u32, u32) {
    map.force_read_attribute::<GeoVertices>(face_id)
        .map(|range| range.0)
        .unwrap_or((0, 0))
}

/// Prepares geo vertices distribution into quadrants for the child cells
fn prepare_geo_vertices_for_subdivision<T: CoordsFloat>(
    map: &CMap2<T>,
    geo_verts: &mut [Vertex2<T>],
    face_darts: &[u32],
    parent_range: (u32, u32),
) -> [(u32, u32); 4] {
    if parent_range.1 == 0 {
        // No vertices to distribute
        return [(0, 0), (0, 0), (0, 0), (0, 0)];
    }

    // Calculate cell bounds for quadrant sorting
    let cell_min = dart_origin(map, face_darts[0]);
    let cell_max = dart_origin(map, face_darts[2]);
    let midpoint = Vertex2::<T>::average(&cell_min, &cell_max);

    distribute_geo_vertices_to_quadrants(geo_verts, parent_range, &midpoint)
}

/// State tracking for subdivision process
struct SubdivisionState {
    going_to_center_prev: u32,
    going_from_center_first: u32,
    dart1_prev: u32,
}

impl SubdivisionState {
    fn new() -> Self {
        Self {
            going_to_center_prev: 0,
            going_from_center_first: 0,
            dart1_prev: 0,
        }
    }
}

/// Performs the actual cell subdivision by splitting edges and creating inner connections
fn perform_cell_subdivision<T: CoordsFloat>(map: &mut CMap2<T>, face_darts: &[u32]) {
    let cell_min = dart_origin(map, face_darts[0]);
    let cell_max = dart_origin(map, face_darts[2]);
    let midpoint = Vertex2::<T>::average(&cell_min, &cell_max);

    let mut state = SubdivisionState::new();

    for (i, &edge_dart) in face_darts.iter().enumerate() {
        // Split boundary edge and create inner edge
        let dart1 = split_boundary_edge(map, edge_dart);
        let (going_to_center, going_from_center) = create_inner_darts(map, edge_dart);

        // Connect edges based on iteration position
        connect_subdivision_edges(
            map,
            i,
            going_to_center,
            going_from_center,
            &mut state,
            dart1,
            &midpoint,
        );
    }
}

/// Connects edges during subdivision based on the current iteration
fn connect_subdivision_edges<T: CoordsFloat>(
    map: &mut CMap2<T>,
    iteration: usize,
    going_to_center: u32,
    going_from_center: u32,
    state: &mut SubdivisionState,
    dart1: u32,
    midpoint: &Vertex2<T>,
) {
    match iteration {
        0 => {
            // First iteration: initialize and set center vertex
            state.going_from_center_first = going_from_center;
            map.force_write_vertex(going_from_center, *midpoint);
        }
        3 => {
            // Last iteration: complete the connections
            map.force_link::<2>(going_from_center, state.going_to_center_prev)
                .unwrap();
            map.force_link::<1>(going_from_center, state.dart1_prev)
                .unwrap();
            map.force_link::<2>(going_to_center, state.going_from_center_first)
                .unwrap();
            map.force_link::<1>(state.going_from_center_first, dart1)
                .unwrap();
        }
        _ => {
            // Middle iterations: connect to previous iteration
            map.force_link::<2>(going_from_center, state.going_to_center_prev)
                .unwrap();
            map.force_link::<1>(going_from_center, state.dart1_prev)
                .unwrap();
        }
    }

    // Update state for next iteration
    state.going_to_center_prev = going_to_center;
    state.dart1_prev = dart1;
}

/// Updates attributes for child cells and returns the face ID that needs further refinement
fn update_child_cell_attributes<T: CoordsFloat>(
    map: &mut CMap2<T>,
    face_darts: &[u32],
    quadrant_ranges: &[(u32, u32); 4],
    parent_face_id: u32,
) -> u32 {
    let parent_refinement_level = map
        .force_read_attribute::<RefinementLevel>(parent_face_id)
        .unwrap()
        .0;

    let mut child_needing_refinement = 0;

    for (i, &dart) in face_darts.iter().enumerate() {
        let child_face_id = map.face_id(dart);

        // Assign geo vertices if this quadrant has any
        if quadrant_ranges[i].1 > 0 {
            map.force_write_attribute::<GeoVertices>(
                child_face_id,
                GeoVertices(quadrant_ranges[i]),
            );
            child_needing_refinement = child_face_id; // One representative child for further refinement
        }

        // Set sibling relationship (circular)
        let next_sibling = map.face_id(face_darts[(i + 1) % 4]);
        map.force_write_attribute::<SiblingDartId>(child_face_id, SiblingDartId(next_sibling));

        // Increment refinement level
        map.force_write_attribute::<RefinementLevel>(
            child_face_id,
            RefinementLevel(parent_refinement_level + 1),
        );
    }

    child_needing_refinement
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
