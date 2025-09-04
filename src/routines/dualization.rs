use crate::routines::helpers::dart_origin;
use honeycomb_core::cmap::{CMap2, CMapBuilder, NULL_DART_ID};
use honeycomb_core::geometry::{CoordsFloat, Vertex2};
use std::collections::{HashMap, HashSet};

/// Iterates on the primal vertices to create a new dual cell around it
/// Uses hashmap for beta2 relations and hashset to avoid duplicating vertices
pub fn dualize_map<T: CoordsFloat>(map: &CMap2<T>) -> CMap2<T> {
    let mut dual_map: CMap2<T> = CMapBuilder::<2, _>::from_n_darts(0).build().unwrap();

    // key : pair of primal vertices, value : dual dart
    let mut dart_hashmap: HashMap<(u32, u32), u32> = HashMap::new();
    // primal_face_id, says if we already created a dual vertex for it
    let mut vertices_hashset: HashSet<u32> = HashSet::new();

    // Step 2: For each vertex in original map, create a dual face
    for vertex_id in map.iter_vertices() {
        // Get all faces around this vertex in order
        let adjacent_primal_vertices = get_adjacent_vertices(map, vertex_id);

        // Skip boundary vertices
        if adjacent_primal_vertices.len() < 3 {
            continue;
        }
        // We should only have quads now
        assert_eq!(adjacent_primal_vertices.len(), 4);

        let new_dart = dual_map.allocate_used_darts(4);

        // If opposite exists we beta2 link to it
        // and write centroid if doesnt exist
        for (i, &remote_primal_vertice) in adjacent_primal_vertices.iter().enumerate() {
            let previous_dart = new_dart + ((i as u32 - 1) % 4); // previous for proper orientation
            let current_new_dart = new_dart + i as u32;

            // Write the centroid
            if !vertices_hashset.contains(&map.face_id(remote_primal_vertice)) {
                let centroid = calculate_face_centroid(map, remote_primal_vertice);

                dual_map.force_write_vertex(previous_dart, centroid);
                vertices_hashset.insert(map.face_id(remote_primal_vertice));
            }

            if dart_hashmap.contains_key(&(map.vertex_id(remote_primal_vertice), vertex_id)) {
                let dual_opposite_dart = *dart_hashmap
                    .get(&(map.vertex_id(remote_primal_vertice), vertex_id))
                    .unwrap();
                dual_map
                    .force_link::<2>(current_new_dart, dual_opposite_dart)
                    .unwrap();
            }

            dart_hashmap.insert(
                (vertex_id, map.vertex_id(remote_primal_vertice)),
                current_new_dart,
            );
        }

        // Connect the darts with beta1 relations to form the dual face
        dual_map.force_link::<1>(new_dart, new_dart + 3).unwrap();
        dual_map
            .force_link::<1>(new_dart + 3, new_dart + 2)
            .unwrap();
        dual_map
            .force_link::<1>(new_dart + 2, new_dart + 1)
            .unwrap();
        dual_map.force_link::<1>(new_dart + 1, new_dart).unwrap();
    }

    dual_map
}

/// Get all faces around a vertex in order
fn get_adjacent_vertices<T: CoordsFloat>(map: &CMap2<T>, vertex_id: u32) -> Vec<u32> {
    assert_ne!(vertex_id, NULL_DART_ID);
    let mut adj_vertices = Vec::new();
    let mut adj_vertices_backward = Vec::new();
    let mut dart = vertex_id;

    adj_vertices.push(map.beta::<1>(dart));

    let mut go_backward = false;
    loop {
        dart = map.beta::<2>(dart);
        dart = map.beta::<1>(dart);

        if dart == NULL_DART_ID {
            go_backward = true;
            break;
        }
        if dart == vertex_id {
            break;
        }

        adj_vertices.push(map.beta::<1>(dart));
    }
    if go_backward {
        dart = vertex_id;
        loop {
            dart = map.beta::<0>(dart);
            dart = map.beta::<2>(dart);

            if dart == NULL_DART_ID {
                break;
            }

            adj_vertices_backward.push(map.beta::<1>(dart));
        }
    }

    adj_vertices_backward.reverse();
    adj_vertices.extend(adj_vertices_backward);

    let adj_vertices_set: std::collections::HashSet<u32> = adj_vertices.iter().cloned().collect();
    assert_eq!(adj_vertices.len(), adj_vertices_set.len());
    assert!(
        !adj_vertices.is_empty(),
        "Vertex {} has less than 1 adj_vertices around it",
        vertex_id
    );

    adj_vertices
}

/// Calculates the centroid of a face by averaging its vertex positions
fn calculate_face_centroid<T: CoordsFloat>(map: &CMap2<T>, vertex_id: u32) -> Vertex2<T> {
    let face_darts: Vec<u32> = map.i_cell::<2>(vertex_id).collect(); //i_cell2 is beta1 orbit

    if face_darts.is_empty() {
        return Vertex2(T::zero(), T::zero());
    }

    let mut sum_x = T::zero();
    let mut sum_y = T::zero();
    let mut count = 0;

    for dart in face_darts {
        let vertex = dart_origin(map, dart);
        sum_x += vertex.x();
        sum_y += vertex.y();
        count += 1;
    }

    if count > 0 {
        let count_t = T::from(count).unwrap();
        Vertex2(sum_x / count_t, sum_y / count_t)
    } else {
        Vertex2(T::zero(), T::zero())
    }
}
