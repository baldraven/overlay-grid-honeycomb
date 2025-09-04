use crate::model::IsIrregular;
use honeycomb_core::{
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
    face_darts
}
