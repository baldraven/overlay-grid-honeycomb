use crate::model::{GeoVertices, IsIrregular, RefinementLevel, SiblingDartId};
use crate::routines::helpers::*;
use honeycomb::prelude::OrbitPolicy;
use honeycomb::{
    core::{
        cmap::{CMap2, NULL_DART_ID},
        geometry::{CoordsFloat, Vertex2},
    },
    stm::{atomically, StmError, Transaction},
};

/// Macro to handle unlink/link operations that may fail, converting errors to StmError::Failure
macro_rules! try_operation {
    ($operation:expr) => {
        match $operation {
            Ok(_) => (),
            Err(_) => return Err(StmError::Failure),
        }
    };
}

/// Result struct for refinement operations
pub struct RefinementResult<T: CoordsFloat> {
    pub children: Vec<u32>,
    pub local_geo_verts: Vec<Vertex2<T>>,
    pub balance_pile: Vec<u32>,
    pub start_end: (usize, usize),
}

/// Fetches the 4 siblings dart to call refine_cell on the 4 of them
pub fn refine_with_pairing<T: CoordsFloat>(
    map: &CMap2<T>,
    working_dart: u32,
    start_dart: u32,
    geo_verts: &[Vertex2<T>],
) -> RefinementResult<T> {
    let mut local_geo_verts = Vec::new();
    let mut children = Vec::new();
    let mut start_end = (0, 0);

    let face1 = map.face_id(working_dart);
    let face2 = map.force_read_attribute::<SiblingDartId>(face1).unwrap().0;
    let face3 = map.force_read_attribute::<SiblingDartId>(face2).unwrap().0;
    let face4 = map.force_read_attribute::<SiblingDartId>(face3).unwrap().0;

    // Mark the future unbalance
    let local_balance_pile = update_balance_pile_for_neighbors(face1, face2, face3, face4, map);

    let result1 = atomically(|trans| refine_cell_tx(trans, map, face1, geo_verts, start_dart));
    let result2 = atomically(|trans| refine_cell_tx(trans, map, face2, geo_verts, start_dart + 16));
    let result3 = atomically(|trans| refine_cell_tx(trans, map, face3, geo_verts, start_dart + 32));
    let result4 = atomically(|trans| refine_cell_tx(trans, map, face4, geo_verts, start_dart + 48));

    let mut results = vec![result1, result2, result3, result4];

    let balance_pile = local_balance_pile;

    //sort the results so that they are in the order of the start_end ranges
    //let mut results = vec![result1, result2, result3, result4];
    results.sort_by_key(|r| r.start_end.0);

    // Combine all
    results.iter().for_each(|result| {
        children.extend(result.children.iter());
        if result.start_end == (0, 0) {
            return;
        } // skip empty results
        local_geo_verts.extend(result.local_geo_verts.iter());
    });

    // Find the overall start and end range from all non-empty results
    let non_empty_ranges: Vec<_> = results
        .iter()
        .filter(|res| res.start_end != (0, 0))
        .map(|res| res.start_end)
        .collect();

    if !non_empty_ranges.is_empty() {
        start_end.0 = non_empty_ranges.first().unwrap().0;
        start_end.1 = non_empty_ranges.last().unwrap().1;
    }

    assert_eq!(
        start_end.1 - start_end.0,
        local_geo_verts.len(),
        "E: local_geo_verts length does not match start_end range"
    );

    RefinementResult {
        children,
        local_geo_verts,
        balance_pile,
        start_end,
    }
}

/// Refines a single cell and returns the child cells that need further refinement (transactional version)
/// to understand better what needs to be done, see README.md
pub fn refine_cell_tx<T: CoordsFloat>(
    trans: &mut Transaction,
    map: &CMap2<T>,
    working_dart: u32,
    geo_verts: &[Vertex2<T>],
    start_dart: u32,
) -> Result<RefinementResult<T>, StmError> {
    let mut local_geo_verts = Vec::new();
    let mut start_end = (0, 0);

    let face_id = working_dart;
    let face_darts = collect_face_darts_tx(trans, map, working_dart)?;
    assert_eq!(face_darts.len(), 4, "E: face darts count is not 4");

    // Step 2: Prepare geo vertices distribution
    let parent_geo_range = map
        .read_attribute::<GeoVertices>(trans, face_id)?
        .unwrap_or(GeoVertices((0, 0)))
        .0;
    map.remove_attribute::<GeoVertices>(trans, face_id)?;

    let mut quadrant_ranges = [(0u32, 0u32); 4];
    if parent_geo_range.1 != 0 {
        let start = parent_geo_range.0 as usize;
        let end = start + parent_geo_range.1 as usize;
        start_end = (start, end);

        local_geo_verts.extend_from_slice(&geo_verts[start..end]);

        // Calculate cell bounds for quadrant sorting
        let cell_min = dart_origin_tx(trans, map, face_darts[0])?;
        let cell_max = dart_origin_tx(trans, map, face_darts[2])?;
        let midpoint = Vertex2::<T>::average(&cell_min, &cell_max);

        // Sort vertices into quadrants: bottom-left, bottom-right, top-right, top-left
        local_geo_verts.sort_by(|a, b| {
            let quadrant_position1 = get_quadrant(a, &midpoint);
            let quadrant_position2 = get_quadrant(b, &midpoint);
            quadrant_position1.cmp(&quadrant_position2)
        });

        // Count vertices in each quadrant
        let mut quadrant_counts = [0u32; 4];
        for vertex in local_geo_verts.iter() {
            let quadrant = get_quadrant(vertex, &midpoint);
            quadrant_counts[quadrant] += 1;
        }

        // Calculate ranges for each quadrant
        let mut current_start = start as u32;

        for i in 0..4 {
            quadrant_ranges[i] = (current_start, quadrant_counts[i]);
            current_start += quadrant_counts[i];
        }

        debug_assert_eq!(
            quadrant_counts.iter().sum::<u32>(),
            parent_geo_range.1,
            "E: quadrant counts do not sum to parent range length"
        );
    }

    // Step 3: Perform the actual cell subdivision
    perform_cell_subdivision_tx(trans, map, &face_darts, start_dart)?;

    // Step 4: Update child cell attributes and find cells needing further refinement
    let children =
        update_child_cell_attributes_tx(trans, map, &face_darts, &quadrant_ranges, face_id)?;

    Ok(RefinementResult {
        children: vec![children],
        local_geo_verts,
        balance_pile: Vec::new(),
        start_end,
    })
}

/// Splits a boundary edge by inserting a new vertex at the midpoint (transactional version)
/// see README.md for visualisation
fn split_boundary_edge_tx<T: CoordsFloat>(
    trans: &mut Transaction,
    map: &CMap2<T>,
    edge_dart: u32,
    start_dart: u32,
) -> Result<u32, StmError> {
    let (n_irregular, next) = canonical_beta1_tx(trans, map, edge_dart)?;
    // If we have an edge that has already been refined because of a refinement of its opposite
    // we need to re-use the extra darts and connect accordingly
    let dart1 = match n_irregular {
        4 => split_boundary_edge_n_4_tx(trans, map, edge_dart)?,
        2 => split_boundary_edge_n_2_tx(trans, map, edge_dart)?,
        1 => split_boundary_edge_n_1_tx(trans, map, edge_dart, next, start_dart)?,
        _ => panic!("E: weirdly irregular ; a cell has probably been over-refined"),
    };
    assert_ne!(dart1, NULL_DART_ID, "E: dart1 is NULL_DART_ID");

    Ok(dart1)
}

/// Simple transactional wrappers for the split_boundary_edge functions
fn split_boundary_edge_n_4_tx<T: CoordsFloat>(
    trans: &mut Transaction,
    map: &CMap2<T>,
    edge_dart: u32,
) -> Result<u32, StmError> {
    let sub_dart1 = map.beta_tx::<1>(trans, edge_dart)?;
    let dart1 = map.beta_tx::<1>(trans, sub_dart1)?;

    // we'll use the pre-existing dart1, so it's not irregular anymore
    map.remove_attribute::<IsIrregular>(trans, dart1)?;
    try_operation!(map.unlink::<1>(trans, sub_dart1));

    Ok(dart1)
}

fn split_boundary_edge_n_2_tx<T: CoordsFloat>(
    trans: &mut Transaction,
    map: &CMap2<T>,
    edge_dart: u32,
) -> Result<u32, StmError> {
    let dart1 = map.beta_tx::<1>(trans, edge_dart)?;

    // we'll use the pre-existing dart1, so it's not irregular anymore
    map.remove_attribute::<IsIrregular>(trans, dart1)?;
    try_operation!(map.unlink::<1>(trans, edge_dart));

    Ok(dart1)
}

fn split_boundary_edge_n_1_tx<T: CoordsFloat>(
    trans: &mut Transaction,
    map: &CMap2<T>,
    edge_dart: u32,
    next: u32,
    start_dart: u32,
) -> Result<u32, StmError> {
    // If we have an opposite we have to also split it and make its new connection, this is how we handle T junctions
    let opposite = map.beta_tx::<2>(trans, edge_dart)?;
    let dart1 = if opposite != NULL_DART_ID {
        map.claim_dart_tx(trans, start_dart)?;
        map.claim_dart_tx(trans, start_dart + 1)?;

        let dart1 = start_dart;
        let dart2 = dart1 + 1;
        map.write_attribute::<IsIrregular>(trans, dart2, IsIrregular(true))?;
        let opposite_next = map.beta_tx::<1>(trans, opposite)?;

        try_operation!(map.unlink::<2>(trans, edge_dart));
        try_operation!(map.unlink::<1>(trans, opposite));
        try_operation!(map.link::<2>(trans, dart1, opposite));
        try_operation!(map.link::<2>(trans, edge_dart, dart2));
        try_operation!(map.link::<1>(trans, opposite, dart2));
        try_operation!(map.link::<1>(trans, dart2, opposite_next));

        dart1
    } else {
        map.claim_dart_tx(trans, start_dart)?;
        start_dart
    };

    try_operation!(map.unlink::<1>(trans, edge_dart));
    try_operation!(map.link::<1>(trans, dart1, next));

    let subdivde_point = Vertex2::<T>::average(
        &dart_origin_tx(trans, map, edge_dart)?,
        &dart_origin_tx(trans, map, next)?,
    );

    map.write_vertex(trans, dart1, subdivde_point)?;

    Ok(dart1)
}

/// Creates inner darts from boundary to center (transactional version)
fn create_inner_darts_tx<T: CoordsFloat>(
    trans: &mut Transaction,
    map: &CMap2<T>,
    edge_dart: u32,
    start_dart: u32,
) -> Result<(u32, u32), StmError> {
    map.claim_dart_tx(trans, start_dart)?;
    map.claim_dart_tx(trans, start_dart + 1)?;

    let going_to_center = start_dart;
    let going_from_center = going_to_center + 1;

    let edge_dart_next = map.beta_tx::<1>(trans, edge_dart)?;

    if edge_dart_next != NULL_DART_ID && !is_regular_tx(trans, map, edge_dart_next)? {
        try_operation!(map.link::<1>(trans, edge_dart_next, going_to_center));
    } else {
        try_operation!(map.link::<1>(trans, edge_dart, going_to_center));
    }
    try_operation!(map.link::<1>(trans, going_to_center, going_from_center));

    Ok((going_to_center, going_from_center))
}

/// Connects edges during subdivision (transactional version)
fn connect_subdivision_edges_tx<T: CoordsFloat>(
    trans: &mut Transaction,
    map: &CMap2<T>,
    iteration: usize,
    going_to_center: u32,
    going_from_center: u32,
    state: &mut SubdivisionState,
    dart1: u32,
    midpoint: &Vertex2<T>,
) -> Result<(), StmError> {
    match iteration {
        0 => {
            // First iteration: initialize and set center vertex
            state.going_from_center_first = going_from_center;
            // let going_from_center_vid = map.vertex_id_tx(trans, going_from_center)?;
            map.write_vertex(trans, going_from_center, *midpoint)?;
        }
        3 => {
            // Last iteration: complete the connections
            try_operation!(map.link::<2>(trans, going_from_center, state.going_to_center_prev));
            try_operation!(map.link::<1>(trans, going_from_center, state.dart1_prev));
            try_operation!(map.link::<2>(trans, going_to_center, state.going_from_center_first));
            try_operation!(map.link::<1>(trans, state.going_from_center_first, dart1));
        }
        _ => {
            // Middle iterations: connect to previous iteration
            try_operation!(map.link::<2>(trans, going_from_center, state.going_to_center_prev));
            try_operation!(map.link::<1>(trans, going_from_center, state.dart1_prev));
        }
    }

    // Update state for next iteration
    state.going_to_center_prev = going_to_center;
    state.dart1_prev = dart1;

    Ok(())
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

/// Performs the actual cell subdivision by splitting edges and creating inner connections (transactional version)
fn perform_cell_subdivision_tx<T: CoordsFloat>(
    trans: &mut Transaction,
    map: &CMap2<T>,
    face_darts: &[u32],
    start_dart: u32,
) -> Result<(), StmError> {
    let cell_min = dart_origin_tx(trans, map, face_darts[0])?;
    let cell_max = dart_origin_tx(trans, map, face_darts[2])?;
    let midpoint = Vertex2::<T>::average(&cell_min, &cell_max);

    let mut state = SubdivisionState::new();

    for (i, &edge_dart) in face_darts.iter().enumerate() {
        assert!(i < 4, "E: face_darts has more than 4 darts");
        // Split boundary edge and create inner edge
        let dart1 = split_boundary_edge_tx(trans, map, edge_dart, start_dart + (i as u32 * 4))?;
        let (going_to_center, going_from_center) =
            create_inner_darts_tx(trans, map, edge_dart, start_dart + (i as u32 * 4) + 2)?;

        // Connect edges based on iteration position
        connect_subdivision_edges_tx(
            trans,
            map,
            i,
            going_to_center,
            going_from_center,
            &mut state,
            dart1,
            &midpoint,
        )?;
    }

    Ok(())
}

/// Updates attributes for child cells and returns the face ID that needs further refinement (transactional version)
fn update_child_cell_attributes_tx<T: CoordsFloat>(
    trans: &mut Transaction,
    map: &CMap2<T>,
    face_darts: &[u32],
    quadrant_ranges: &[(u32, u32); 4],
    parent_face_id: u32,
) -> Result<u32, StmError> {
    let parent_refinement_level = map
        .read_attribute::<RefinementLevel>(trans, parent_face_id)?
        .unwrap()
        .0;

    let mut child_needing_refinement = 0;

    for (i, &dart) in face_darts.iter().enumerate() {
        let child_face_id = map.face_id_tx(trans, dart)?;

        // Assign geo vertices if this quadrant has any
        if quadrant_ranges[i].1 > 0 {
            map.write_attribute::<GeoVertices>(
                trans,
                child_face_id,
                GeoVertices(quadrant_ranges[i]),
            )?;
            child_needing_refinement = child_face_id; // One representative child for further refinement
        }

        // Set sibling relationship (circular)
        let next_sibling = map.face_id_tx(trans, face_darts[(i + 1) % 4])?;
        map.write_attribute::<SiblingDartId>(trans, child_face_id, SiblingDartId(next_sibling))?;

        // Increment refinement level
        map.write_attribute::<RefinementLevel>(
            trans,
            child_face_id,
            RefinementLevel(parent_refinement_level + 1),
        )?;
    }

    Ok(child_needing_refinement)
}

//--------------- NON TX

fn update_balance_pile_for_neighbors<T: CoordsFloat>(
    face1: u32,
    face2: u32,
    face3: u32,
    face4: u32,
    map: &CMap2<T>,
) -> Vec<u32> {
    let current_depth = map
        .force_read_attribute::<RefinementLevel>(face1)
        .unwrap()
        .0;

    let mut balance_stack = Vec::new();

    for i in 0..4 {
        let face = match i {
            0 => face1,
            1 => face2,
            2 => face3,
            3 => face4,
            _ => unreachable!(),
        };
        let orbit = map.orbit(OrbitPolicy::Face, face).collect::<Vec<u32>>();

        for dart in orbit {
            let opposite = map.beta::<2>(dart);
            if opposite == NULL_DART_ID {
                continue;
            }

            let opposite_face = map.face_id(opposite);
            let opposite_depth = map
                .force_read_attribute::<RefinementLevel>(opposite_face)
                .unwrap()
                .0;

            if current_depth > opposite_depth {
                balance_stack.push(opposite_face);
            }
        }
    }

    balance_stack.sort();
    balance_stack.dedup();

    balance_stack
}

/// Performs the actual cell subdivision by splitting edges and creating inner connections
fn perform_cell_subdivision<T: CoordsFloat>(map: &CMap2<T>, face_darts: &[u32], start_dart: u32) {
    let cell_min = dart_origin(map, face_darts[0]);
    let cell_max = dart_origin(map, face_darts[2]);
    let midpoint = Vertex2::<T>::average(&cell_min, &cell_max);

    let mut state = SubdivisionState::new();

    for (i, &edge_dart) in face_darts.iter().enumerate() {
        // Split boundary edge and create inner edge
        let dart1 = split_boundary_edge(map, edge_dart, start_dart + (i as u32 * 4));
        let (going_to_center, going_from_center) =
            create_inner_darts(map, edge_dart, start_dart + (i as u32 * 4) + 2);

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
    map: &CMap2<T>,
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

/// Refines a single cell and returns the child cells that need further refinement
/// to understand better what needs to be done, see README.md
pub fn refine_cell<T: CoordsFloat>(
    map: &CMap2<T>,
    working_dart: u32,
    geo_verts: &[Vertex2<T>],
    start_dart: u32,
) -> RefinementResult<T> {
    let mut local_geo_verts = Vec::new();
    let mut start_end = (0, 0);

    let face_id = map.face_id(working_dart);
    let face_darts = collect_face_darts(map, working_dart);
    assert_eq!(face_darts.len(), 4, "E: face darts count is not 4");

    // Step 2: Prepare geo vertices distribution
    let parent_geo_range = map
        .force_read_attribute::<GeoVertices>(face_id)
        .unwrap_or(GeoVertices((0, 0)))
        .0;
    map.force_remove_attribute::<GeoVertices>(face_id);

    let mut quadrant_ranges = [(0u32, 0u32); 4];
    if parent_geo_range.1 != 0 {
        let start = parent_geo_range.0 as usize;
        let end = start + parent_geo_range.1 as usize;
        start_end = (start, end);

        local_geo_verts.clear();
        local_geo_verts.extend_from_slice(&geo_verts[start..end]);

        // Calculate cell bounds for quadrant sorting
        let cell_min = dart_origin(map, face_darts[0]);
        let cell_max = dart_origin(map, face_darts[2]);
        let midpoint = Vertex2::<T>::average(&cell_min, &cell_max);

        // Sort vertices into quadrants: bottom-left, bottom-right, top-right, top-left
        local_geo_verts.sort_by(|a, b| {
            let quadrant_position1 = get_quadrant(a, &midpoint);
            let quadrant_position2 = get_quadrant(b, &midpoint);
            quadrant_position1.cmp(&quadrant_position2)
        });

        // Count vertices in each quadrant
        let mut quadrant_counts = [0u32; 4];
        for vertex in local_geo_verts.iter() {
            let quadrant = get_quadrant(vertex, &midpoint);
            quadrant_counts[quadrant] += 1;
        }

        // Calculate ranges for each quadrant
        let mut current_start = start as u32;

        for i in 0..4 {
            quadrant_ranges[i] = (current_start, quadrant_counts[i]);
            current_start += quadrant_counts[i];
        }

        assert_eq!(
            quadrant_counts.iter().sum::<u32>(),
            parent_geo_range.1,
            "E: quadrant counts do not sum to parent range length"
        );
    }

    // Step 3: Perform the actual cell subdivision
    perform_cell_subdivision(map, &face_darts, start_dart);

    // Step 4: Update child cell attributes and find cells needing further refinement
    let children = update_child_cell_attributes(map, &face_darts, &quadrant_ranges, face_id);

    RefinementResult {
        children: vec![children],
        local_geo_verts,
        balance_pile: Vec::new(),
        start_end,
    }
}

/// Updates attributes for child cells and returns the face ID that needs further refinement
fn update_child_cell_attributes<T: CoordsFloat>(
    map: &CMap2<T>,
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

fn split_boundary_edge_n_1<T: CoordsFloat>(
    map: &CMap2<T>,
    edge_dart: u32,
    next: u32,
    start_dart: u32,
) -> u32 {
    // If we have an opposite we have to also split it and make its new connection, this is how we handle T junctions
    let opposite = map.beta::<2>(edge_dart);
    let dart1 = if opposite != NULL_DART_ID {
        atomically(|trans| {
            map.claim_dart_tx(trans, start_dart)?;
            map.claim_dart_tx(trans, start_dart + 1)?;
            Ok(())
        });
        let dart1 = start_dart;
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
        atomically(|trans| {
            map.claim_dart_tx(trans, start_dart)?;
            Ok(())
        });
        start_dart
    };

    map.force_unlink::<1>(edge_dart).unwrap();
    map.force_link::<1>(dart1, next).unwrap();

    let subdivde_point =
        Vertex2::<T>::average(&dart_origin(map, edge_dart), &dart_origin(map, next));
    map.force_write_vertex(dart1, subdivde_point);

    dart1
}

fn split_boundary_edge_n_4<T: CoordsFloat>(map: &CMap2<T>, edge_dart: u32) -> u32 {
    let sub_dart1 = map.beta::<1>(edge_dart);
    let dart1 = map.beta::<1>(sub_dart1);

    // we'll use the pre-existing dart1, so it's not irregular anymore
    map.force_remove_attribute::<IsIrregular>(dart1);
    map.force_unlink::<1>(sub_dart1).unwrap();

    dart1
}

fn split_boundary_edge_n_2<T: CoordsFloat>(map: &CMap2<T>, edge_dart: u32) -> u32 {
    let dart1 = map.beta::<1>(edge_dart);

    // we'll use the pre-existing dart1, so it's not irregular anymore
    map.force_remove_attribute::<IsIrregular>(dart1);
    map.force_unlink::<1>(edge_dart).unwrap();

    dart1
}

/// Creates inner darts from boundary to center
fn create_inner_darts<T: CoordsFloat>(
    map: &CMap2<T>,
    edge_dart: u32,
    start_dart: u32,
) -> (u32, u32) {
    atomically(|trans| {
        map.claim_dart_tx(trans, start_dart)?;
        map.claim_dart_tx(trans, start_dart + 1)?;
        Ok(())
    });
    let going_to_center = start_dart;
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

/// Splits a boundary edge by inserting a new vertex at the midpoint
/// see README.md for visualisation
fn split_boundary_edge<T: CoordsFloat>(map: &CMap2<T>, edge_dart: u32, start_dart: u32) -> u32 {
    let (n_irregular, next) = canonical_beta1(map, edge_dart);
    // If we have an edge that has already been refined because of a refinement of its opposite
    // we need to re-use the extra darts and connect accordingly
    let dart1 = match n_irregular {
        4 => split_boundary_edge_n_4(map, edge_dart),
        2 => split_boundary_edge_n_2(map, edge_dart),
        1 => split_boundary_edge_n_1(map, edge_dart, next, start_dart),
        _ => panic!("E: weirdly irregular ; a cell has probably been over-refined"),
    };
    assert_ne!(dart1, NULL_DART_ID, "E: dart1 is NULL_DART_ID");

    dart1
}
