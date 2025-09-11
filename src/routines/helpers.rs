use crate::model::IsIrregular;
use honeycomb::core::{
    cmap::{CMap2, NULL_DART_ID},
    geometry::{CoordsFloat, Vertex2},
};

pub fn dart_origin<T: CoordsFloat>(map: &CMap2<T>, dart: u32) -> Vertex2<T> {
    assert_ne!(dart, NULL_DART_ID);

    map.force_read_vertex(map.vertex_id(dart)).unwrap()
}

pub fn canonical_beta1<T: CoordsFloat>(map: &CMap2<T>, dart: u32) -> (u32, u32) {
    assert_ne!(dart, NULL_DART_ID);

    let mut steps = 1;
    let mut current_dart = map.beta::<1>(dart);

    while !is_regular(map, current_dart) {
        steps += 1;
        current_dart = map.beta::<1>(current_dart);
    }

    (steps, current_dart)
}

pub fn is_regular<T: CoordsFloat>(map: &CMap2<T>, dart: u32) -> bool {
    assert_ne!(dart, NULL_DART_ID);

    let attr = map.force_read_attribute::<IsIrregular>(dart);
    attr.is_none()
}

/// Collects all darts non-irregular in a face by following the beta1 relation
/// Make sure geometrically that first dart is south-west, useful for the geo_vertices distribution
pub fn collect_face_darts<T: CoordsFloat>(map: &CMap2<T>, start_dart: u32) -> Vec<u32> {
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

    // order the darts so the first face_dart is in SW quadrant
    let pos0 = dart_origin(map, face_darts[0]);
    let pos2 = dart_origin(map, face_darts[2]);
    let midpoint = Vertex2::<T>::average(&pos0, &pos2);
    let mut face_darts_ordered = vec![0u32; 4];
    for dart in face_darts.iter() {
        let pos = dart_origin(map, *dart);
        let quadrant = get_quadrant(&pos, &midpoint);
        face_darts_ordered[quadrant] = *dart;
    }

    face_darts_ordered
}

/// Determines which quadrant a vertex belongs to relative to a midpoint
/// Returns: 0=bottom-left, 1=bottom-right, 2=top-right, 3=top-left
pub fn get_quadrant<T: CoordsFloat>(vertex: &Vertex2<T>, midpoint: &Vertex2<T>) -> usize {
    let x_right = vertex.x() >= midpoint.x();
    let y_top = vertex.y() >= midpoint.y();

    match (x_right, y_top) {
        (false, false) => 0, // bottom-left
        (true, false) => 1,  // bottom-right
        (true, true) => 2,   // top-right
        (false, true) => 3,  // top-left
    }
}
