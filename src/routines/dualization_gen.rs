use honeycomb_core::geometry::{CoordsFloat, Vertex2};
use honeycomb_core::cmap::{CMap2, CMapBuilder};
use std::collections::HashMap;
use crate::routines::helpers::dart_origin;

pub fn dualize_map<T: CoordsFloat>(map: &CMap2<T>) -> CMap2<T> {
    // Create a new CMap2 for the dual
    let mut dual_map: CMap2<T> = CMapBuilder::<2, _>::from_n_darts(0)
        .build()
        .unwrap();
    
    // Step 1: Create vertices in dual map from face centroids of original map
    let face_to_vertex_map = create_dual_vertices(&mut dual_map, map);
    
    // Step 2: Create faces in dual map from vertices of original map
    create_dual_faces(&mut dual_map, map, &face_to_vertex_map);
    
    dual_map
}

/// Creates vertices in the dual map from face centroids of the original map
fn create_dual_vertices<T: CoordsFloat>(
    dual_map: &mut CMap2<T>,
    original_map: &CMap2<T>,
) -> HashMap<u32, u32> {
    let mut face_to_vertex_map = HashMap::new();
    
    for face_id in original_map.iter_faces() {
        // Calculate the centroid of this face
        let centroid = calculate_face_centroid(original_map, face_id);
        
        // Create a new vertex in the dual map at the centroid position
        let dual_vertex_dart = dual_map.allocate_used_darts(1);
        dual_map.force_write_vertex(dual_vertex_dart, centroid);
        
        // Map original face ID to dual vertex dart
        face_to_vertex_map.insert(face_id, dual_vertex_dart);
    }
    
    face_to_vertex_map
}

/// Calculates the centroid of a face by averaging its vertex positions
fn calculate_face_centroid<T: CoordsFloat>(map: &CMap2<T>, face_id: u32) -> Vertex2<T> {
    let face_darts: Vec<u32> = map.i_cell::<1>(face_id).collect();
    
    if face_darts.is_empty() {
        return Vertex2(T::zero(), T::zero());
    }
    
    let mut sum_x = T::zero();
    let mut sum_y = T::zero();
    let mut count = 0;
    
    for dart in face_darts {
        let vertex = dart_origin(map, dart);
        sum_x = sum_x + vertex.x();
        sum_y = sum_y + vertex.y();
        count += 1;
    }
    
    if count > 0 {
        let count_t = T::from(count).unwrap();
        Vertex2(sum_x / count_t, sum_y / count_t)
    } else {
        Vertex2(T::zero(), T::zero())
    }
}

/// Creates faces in the dual map from vertices of the original map
fn create_dual_faces<T: CoordsFloat>(
    dual_map: &mut CMap2<T>,
    original_map: &CMap2<T>,
    face_to_vertex_map: &HashMap<u32, u32>,
) {
    // Collect all vertex representatives to avoid processing the same vertex multiple times
    let mut processed_vertices = std::collections::HashSet::new();
    
    // For each vertex in the original map, create a face in the dual map
    for dart in original_map.iter_faces().flat_map(|face| original_map.i_cell::<1>(face)) {
        let vertex_id = original_map.vertex_id(dart);
        
        if processed_vertices.contains(&vertex_id) {
            continue;
        }
        processed_vertices.insert(vertex_id);
        
        // Get all faces around this vertex
        let surrounding_faces = get_faces_around_vertex(original_map, dart);
        
        if surrounding_faces.len() >= 3 {
            // Create a face in dual map connecting the corresponding dual vertices
            create_dual_face_from_surrounding_faces(dual_map, &surrounding_faces, face_to_vertex_map);
        }
    }
}

/// Gets all faces surrounding a vertex
fn get_faces_around_vertex<T: CoordsFloat>(map: &CMap2<T>, vertex_dart: u32) -> Vec<u32> {
    let mut faces = Vec::new();
    let vertex_darts: Vec<u32> = map.i_cell::<0>(vertex_dart).collect();
    
    for dart in vertex_darts {
        let face_id = map.face_id(dart);
        if !faces.contains(&face_id) {
            faces.push(face_id);
        }
    }
    
    faces
}

/// Creates a dual face by connecting dual vertices corresponding to surrounding faces
fn create_dual_face_from_surrounding_faces<T: CoordsFloat>(
    dual_map: &mut CMap2<T>,
    surrounding_faces: &[u32],
    face_to_vertex_map: &HashMap<u32, u32>,
) {
    if surrounding_faces.len() < 3 {
        return;
    }
    
    // Get the dual vertex darts for each surrounding face
    let mut dual_vertex_positions = Vec::new();
    for &face_id in surrounding_faces {
        if let Some(&dual_vertex_dart) = face_to_vertex_map.get(&face_id) {
            let pos = dart_origin(dual_map, dual_vertex_dart);
            dual_vertex_positions.push(pos);
        }
    }
    
    if dual_vertex_positions.len() < 3 {
        return;
    }
    
    // Create edges connecting these dual vertices to form a face
    create_face_from_vertex_positions(dual_map, &dual_vertex_positions);
}

/// Creates a face in the dual map by connecting vertices at the given positions
fn create_face_from_vertex_positions<T: CoordsFloat>(dual_map: &mut CMap2<T>, vertex_positions: &[Vertex2<T>]) {
    if vertex_positions.len() < 3 {
        return;
    }
    
    let n = vertex_positions.len();
    let mut edge_darts = Vec::new();
    
    // Create edge darts for the face
    for i in 0..n {
        let dart1 = dual_map.allocate_used_darts(2);
        let dart2 = dart1 + 1;
        edge_darts.push((dart1, dart2));
        
        // Set the vertex positions for this edge
        let current_pos = vertex_positions[i];
        let next_pos = vertex_positions[(i + 1) % n];
        
        dual_map.force_write_vertex(dart1, current_pos);
        dual_map.force_write_vertex(dart2, next_pos);
    }
    
    // Link the edges to form a face
    for i in 0..n {
        let (current_dart1, current_dart2) = edge_darts[i];
        let (next_dart1, _) = edge_darts[(i + 1) % n];
        let (_, prev_dart2) = edge_darts[(i + n - 1) % n];
        
        // Beta-1 connections (next edge in face)
        dual_map.force_link::<1>(current_dart1, next_dart1).unwrap();
        dual_map.force_link::<1>(prev_dart2, current_dart2).unwrap();
        
        // Beta-2 connections (opposite edge)
        dual_map.force_link::<2>(current_dart1, current_dart2).unwrap();
    }
}