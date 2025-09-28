use honeycomb::core::geometry::CoordsFloat;

use crate::routines::helpers::dart_origin;
use honeycomb::core::cmap::{CMap2, NULL_DART_ID};
use honeycomb::core::geometry::Vertex2;

pub fn _print_map_info<T: CoordsFloat>(map: &CMap2<T>) {
    let num_vertices = map.iter_vertices().count() as u32;
    let num_edges = map.iter_edges().count() as u32;
    let num_faces = map.iter_faces().count() as u32;

    println!("Map Information:");
    println!("Number of Vertices: {}", num_vertices);
    println!("Number of Edges: {}", num_edges);
    println!("Number of Faces: {}", num_faces);

    // Print all vertices
    for vertex_id in map.iter_vertices() {
        let vertex = map.force_read_vertex(map.vertex_id(vertex_id));
        if vertex.is_some() {
            println!("Vertex {:?}", vertex);
        } else {
            println!("Vertex with ID {} has no associated geometry.", vertex_id);
        }
    }
}

/// Debugging function to create a triangle around a given dart in 2D
/// Creates a small triangle with the dart's origin as one vertex and the beta-1 dart's origin as another
pub fn _triangle_around_dart<T: CoordsFloat>(map: &mut CMap2<T>, dart_id: u32) {
    assert_ne!(dart_id, NULL_DART_ID);
    let base_vertex = dart_origin(map, dart_id);
    let next_vertex = dart_origin(map, map.beta::<1>(dart_id));

    // Allocate 3 darts for the triangle
    let triangle_dart = map.allocate_used_darts(3);

    // Compute the direction and length of the dart
    let dir = Vertex2(next_vertex.0 - base_vertex.0, next_vertex.1 - base_vertex.1);
    let len = (dir.0 * dir.0 + dir.1 * dir.1).sqrt();
    let frac = T::from(0.15).unwrap(); // Use 15% of the dart length for offset

    // Create the third vertex perpendicular to the dart direction
    // Rotate the direction vector 90 degrees counterclockwise
    let perp = Vertex2(-dir.1, dir.0); // 90° rotation
    let perp_normalized = if len > T::from(0.0).unwrap() {
        Vertex2(perp.0 / len, perp.1 / len)
    } else {
        Vertex2(T::from(0.0).unwrap(), T::from(1.0).unwrap()) // fallback
    };

    let third_vertex = Vertex2(
        base_vertex.0 + frac * len * perp_normalized.0,
        base_vertex.1 + frac * len * perp_normalized.1,
    );

    // Set vertices for the triangle
    map.force_write_vertex(triangle_dart, base_vertex);
    map.force_write_vertex(triangle_dart + 1, next_vertex);
    map.force_write_vertex(triangle_dart + 2, third_vertex);

    // Link the triangle (beta-1 connections)
    map.force_link::<1>(triangle_dart, triangle_dart + 1)
        .unwrap();
    map.force_link::<1>(triangle_dart + 1, triangle_dart + 2)
        .unwrap();
    map.force_link::<1>(triangle_dart + 2, triangle_dart)
        .unwrap();
}
