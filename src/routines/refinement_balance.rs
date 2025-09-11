use std::collections::{HashSet, VecDeque};

use crate::{
    model::{RefinementLevel, SiblingDartId},
    routines::{edge_splitting::*, pre_processing::*},
};
use honeycomb::core::{
    cmap::{CMap2, NULL_DART_ID},
    geometry::{CoordsFloat, Vertex2},
};

/// Performs refinement of cells that contain geo vertices
/// Processes one level at a time to ensure proper ordering
pub fn refinement<T: CoordsFloat>(
    map: &mut CMap2<T>,
    geo_verts: &mut [Vertex2<T>],
    max_depth: u32,
) -> Result<(), OverlayGridError> {
    let mut balance_pile = VecDeque::new();
    let mut current_level = VecDeque::new();
    let to_refine = refine_cell(map, 1, geo_verts);
    current_level.push_back(to_refine);

    // Process each refinement level
    for depth in 0..max_depth {
        assert_ne!(current_level.len(), 0);

        println!(
            "Processing refinement level {} with {} pack of cells to refine",
            depth,
            current_level.len()
        );

        let mut next_level = VecDeque::new();

        // Refine all cells at current level
        // TODO: parallelize this
        while let Some(working_dart) = current_level.pop_front() {
            let parent_cells = refine_with_pairing(map, working_dart, geo_verts, &mut balance_pile);

            // Add children to next level
            for child_cell in parent_cells {
                if child_cell != NULL_DART_ID {
                    next_level.push_back(child_cell);
                }
            }
        }

        balancing(map, &mut balance_pile, geo_verts);

        current_level = next_level;
    }

    Ok(())
}

fn balancing<T: CoordsFloat>(
    map: &mut CMap2<T>,
    balance_pile: &mut VecDeque<u32>,
    geo_verts: &mut [Vertex2<T>],
) {
    let mut seen_faces: HashSet<u32> = HashSet::new();

    while !balance_pile.is_empty() {
        let face_dart = map.face_id(balance_pile.pop_front().unwrap());
        if seen_faces.contains(&face_dart) {
            continue;
        }

        if ensure_not_balanced(map, face_dart) {
            // gather the 4 sibling faces (first ring) and insert them into seen_faces
            let mut first_ring = Vec::with_capacity(4);
            let mut current_dart = map.face_id(face_dart);
            first_ring.push(current_dart);
            for _ in 1..4 {
                current_dart = map
                    .force_read_attribute::<SiblingDartId>(current_dart)
                    .unwrap()
                    .0;
                first_ring.push(current_dart);
            }

            refine_with_pairing(map, face_dart, geo_verts, balance_pile);

            // gather the next-level siblings (3 successors for each of the 4 faces => 12 faces)
            let mut second_ring = Vec::with_capacity(12);
            for &f in &first_ring {
                let mut current_dart = f;
                for _ in 0..3 {
                    current_dart = map
                        .force_read_attribute::<SiblingDartId>(current_dart)
                        .unwrap()
                        .0;
                    second_ring.push(current_dart);
                }
            }
            for &f in &first_ring {
                seen_faces.insert(f);
            }
            for &f in &second_ring {
                seen_faces.insert(f);
            }
        }
    }
}

fn ensure_not_balanced<T: CoordsFloat>(map: &CMap2<T>, working_dart: u32) -> bool {
    let face_id = map.face_id(working_dart);

    // Get the current cell's refinement level
    let current_level = map
        .force_read_attribute::<RefinementLevel>(face_id)
        .unwrap()
        .0;

    for dart in map.i_cell::<2>(face_id) {
        let opposite = map.beta::<2>(dart);
        if opposite == NULL_DART_ID {
            continue; // Skip boundary edges
        }

        let opposite_face_id = map.face_id(opposite);
        let opposite_level = map
            .force_read_attribute::<RefinementLevel>(opposite_face_id)
            .unwrap()
            .0;

        // If there's a refinement level difference > 1, the cell is unbalanced
        if (current_level as i32 - opposite_level as i32).abs() > 1 {
            return true; // Cell needs balancing
        }
    }

    println!("is already balanced");
    false
}
