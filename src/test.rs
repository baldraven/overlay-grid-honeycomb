use super::*;
use crate::{model::Geometry2, routines::compute_overlapping_grid_size};
use honeycomb_core::{cmap::CMapBuilder, geometry::Vertex2};

#[test]
fn test_collect_face_darts_quad() {
    let map: CMap2<f64> = CMapBuilder::<2, _>::unit_grid(1).build().unwrap();
    let face_darts = collect_face_darts(&map, 1);

    assert_eq!(face_darts.len(), 4);
    assert_eq!(face_darts, vec![1, 2, 3, 4]);
}

#[test]
fn test_collect_face_darts_rectangle() {
    let map: CMap2<f64> = CMapBuilder::<2, _>::unit_grid(1).build().unwrap();
    let face_darts = collect_face_darts(&map, 1);

    assert_eq!(face_darts.len(), 4);
    assert_eq!(face_darts, vec![1, 2, 3, 4]);
}

#[test]
fn test_collect_face_darts_different_start() {
    let map: CMap2<f64> = CMapBuilder::<2, _>::unit_grid(1).build().unwrap();
    let face_darts_from_2 = collect_face_darts(&map, 2);

    assert_eq!(face_darts_from_2.len(), 4);
    assert_eq!(face_darts_from_2, vec![2, 3, 4, 1]);
}

#[test]
fn test_subdivision_state_default() {
    let state = SubdivisionState::default();

    // Default values should be 0
    assert_eq!(state.going_to_center_prev, 0);
    assert_eq!(state.going_from_center_first, 0);
    assert_eq!(state.dart1_prev, 0);
}

#[test]
fn test_subdivision_state_update() {
    let mut state = SubdivisionState::default();

    state.update_for_iteration(10, 20);

    assert_eq!(state.going_to_center_prev, 10);
    assert_eq!(state.dart1_prev, 20);
    // going_from_center_first should remain unchanged
    assert_eq!(state.going_from_center_first, 0);
}

#[test]
fn test_subdivision_state_multiple_updates() {
    let mut state = SubdivisionState::default();

    // First update
    state.update_for_iteration(10, 20);
    assert_eq!(state.going_to_center_prev, 10);
    assert_eq!(state.dart1_prev, 20);

    // Second update should overwrite previous values
    state.update_for_iteration(30, 40);
    assert_eq!(state.going_to_center_prev, 30);
    assert_eq!(state.dart1_prev, 40);
}

#[test]
fn test_subdivision_state_first_iteration_simulation() {
    let mut state = SubdivisionState::default();

    // Simulate first iteration behavior
    state.going_from_center_first = 100;
    state.update_for_iteration(10, 20);

    assert_eq!(state.going_from_center_first, 100); // Should remain set
    assert_eq!(state.going_to_center_prev, 10);
    assert_eq!(state.dart1_prev, 20);
}

#[test]
fn test_error_handling_create_inner_edge_minimal() {
    let mut map: CMap2<f64> = CMapBuilder::<2, _>::from_n_darts(1).build().unwrap();

    // Trying to create inner edge from a dart that's not properly set up
    let result = create_inner_edge(&mut map, 1);
    // This should work (just allocate darts and link), but let's verify the behavior
    assert!(result.is_ok());
}

#[test]
fn test_collect_face_darts_single_dart_loop() {
    let map: CMap2<f64> = CMapBuilder::<2, _>::from_n_darts(1).build().unwrap();
    map.force_link::<1>(1, 1).unwrap(); // Self-loop

    let face_darts = collect_face_darts(&map, 1);

    assert_eq!(face_darts.len(), 1);
    assert_eq!(face_darts, vec![1]);
}

#[test]
fn test_compute_overlapping_grid_size_valid() {
    let geometry: Geometry2<f64> = Geometry2 {
        vertices: vec![
            Vertex2::from((0.0, 0.0)),
            Vertex2::from((2.0, 0.0)),
            Vertex2::from((1.0, 3.0)),
            Vertex2::from((-1.0, -1.0)),
        ],
        segments: vec![(0, 1), (1, 2), (2, 3), (3, 0)],
    };

    let result = compute_overlapping_grid_size(&geometry);
    assert!(result.is_ok());

    let [max_x, max_y, min_x, min_y] = result.unwrap();
    assert_eq!(max_x, 2.0);
    assert_eq!(max_y, 3.0);
    assert_eq!(min_x, -1.0);
    assert_eq!(min_y, -1.0);
}

// Integration test for helper functions working together

#[test]
fn test_helper_functions_integration() {
    let mut map = CMapBuilder::<2, _>::unit_grid(1).build().unwrap();

    // Set vertices for a proper rectangle
    map.force_write_vertex(1, Vertex2::from((0.0, 0.0)));
    map.force_write_vertex(2, Vertex2::from((1.0, 0.0)));
    map.force_write_vertex(3, Vertex2::from((1.0, 1.0)));
    map.force_write_vertex(4, Vertex2::from((0.0, 1.0)));

    // Test the workflow that's used in the main overlay_grid function
    let working_dart = 1;
    let face_darts = collect_face_darts(&map, working_dart);
    assert_eq!(face_darts.len(), 4);

    let mut state = SubdivisionState::default();

    // Simulate the first few iterations of the subdivision process
    for (i, &edge_dart) in face_darts.iter().take(2).enumerate() {
        // Split boundary edge
        let dart1 = split_boundary_edge(&mut map, edge_dart).unwrap();

        // Create inner edge
        let (going_to_center, going_from_center) = create_inner_edge(&mut map, edge_dart).unwrap();

        if i == 0 {
            // First iteration: set the center vertex
            let center = Vertex2::from((0.5, 0.5)); // Center of rectangle
            map.force_write_vertex(going_from_center, center);
            state.going_from_center_first = going_from_center;
        } else {
            // Subsequent iterations: link with previous
            assert!(map
                .force_link::<2>(going_from_center, state.going_to_center_prev)
                .is_ok());
            assert!(map
                .force_link::<1>(going_from_center, state.dart1_prev)
                .is_ok());
        }

        state.update_for_iteration(going_to_center, dart1);
    }

    // Verify that the state was properly updated
    assert_ne!(state.going_to_center_prev, 0);
    assert_ne!(state.dart1_prev, 0);
}
#[test]
fn test_get_quadrant_functionality() {
    let midpoint = Vertex2::<f64>::from((0.5, 0.5));

    // Test all four quadrants
    assert_eq!(get_quadrant(&Vertex2::from((0.2, 0.2)), &midpoint), 0); // bottom-left
    assert_eq!(get_quadrant(&Vertex2::from((0.8, 0.2)), &midpoint), 1); // bottom-right
    assert_eq!(get_quadrant(&Vertex2::from((0.8, 0.8)), &midpoint), 2); // top-right
    assert_eq!(get_quadrant(&Vertex2::from((0.2, 0.8)), &midpoint), 3); // top-left

    // Test edge cases (on the midpoint lines)
    assert_eq!(get_quadrant(&Vertex2::from((0.5, 0.2)), &midpoint), 1); // on vertical midline, bottom
    assert_eq!(get_quadrant(&Vertex2::from((0.2, 0.5)), &midpoint), 3); // on horizontal midline, left
}

#[test]
fn test_distribute_geo_vertices_empty_range() {
    let mut vertices = vec![
        Vertex2::<f64>::from((0.1, 0.1)),
        Vertex2::<f64>::from((0.9, 0.1)),
    ];
    let midpoint = Vertex2::<f64>::from((0.5, 0.5));

    // Test with empty range
    let result = distribute_geo_vertices_to_quadrants(&mut vertices, (0, 0), &midpoint);
    assert_eq!(result, [(0, 0), (0, 0), (0, 0), (0, 0)]);
}

#[test]
fn test_distribute_geo_vertices_single_vertex() {
    let mut vertices = vec![
        Vertex2::<f64>::from((0.2, 0.2)), // bottom-left quadrant
    ];
    let midpoint = Vertex2::<f64>::from((0.5, 0.5));

    let result = distribute_geo_vertices_to_quadrants(&mut vertices, (0, 1), &midpoint);

    // All vertices should be in bottom-left quadrant (index 0)
    assert_eq!(result[0], (0, 1)); // bottom-left: 1 vertex
    assert_eq!(result[1], (1, 0)); // bottom-right: 0 vertices
    assert_eq!(result[2], (1, 0)); // top-right: 0 vertices
    assert_eq!(result[3], (1, 0)); // top-left: 0 vertices
}
