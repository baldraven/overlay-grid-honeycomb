use std::collections::HashSet;

use crate::{
    model::SiblingDartId,
    routines::{edge_splitting::*, pre_processing::*},
};
use honeycomb::core::{
    cmap::{CMap2, NULL_DART_ID},
    geometry::{CoordsFloat, Vertex2},
};

use rayon::prelude::*;

pub fn refinement<T: CoordsFloat>(
    map: &mut CMap2<T>,
    geo_verts: &mut [Vertex2<T>],
    max_depth: u32,
) -> Result<(), OverlayGridError> {
    let mut balance_pile = Vec::new();
    let mut darts_to_refine = Vec::new();

    // Start by refining one cell to initialize the process
    let dart = map.allocate_unused_darts(16);
    let result = refine_cell(map, 1, geo_verts, dart);
    let to_refine = result.children[0];
    darts_to_refine.push(to_refine);
    if !result.local_geo_verts.is_empty() {
        geo_verts[result.start_end.0..result.start_end.1].copy_from_slice(&result.local_geo_verts);
    }

    // Main loop
    for depth in 0..max_depth {
        let nb_allocations = 64 * darts_to_refine.len();
        let new_dart = map.allocate_unused_darts(nb_allocations);
        let chunks = (0..darts_to_refine.len()).map(|i| new_dart + (i as u32 * 64));

        println!(
            "Processing refinement level {} with {} pack of cells to refine",
            depth,
            darts_to_refine.len()
        );

        let current_level_zipped: Vec<(u32, u32)> =
            darts_to_refine.into_iter().zip(chunks).collect();

        let mut next_darts_to_refine = Vec::new();
        parallel_refine_with_pairing(
            map,
            current_level_zipped.into_par_iter(),
            &mut next_darts_to_refine,
            &mut balance_pile,
            geo_verts,
        );

        balancing(map, &mut balance_pile, geo_verts);

        darts_to_refine = next_darts_to_refine;
    }

    Ok(())
}

/// Sanitizes the balance pile to remove duplicates and keep only one representative from each sibling group
fn sanitize_balance_pile<T: CoordsFloat>(map: &CMap2<T>, balance_pile: &Vec<u32>) -> Vec<u32> {
    let mut processed_faces = HashSet::new();
    let mut sanitized = Vec::new();

    for &face_dart in balance_pile {
        if processed_faces.contains(&face_dart) {
            continue;
        }

        // Get all 4 sibling faces in this group
        let mut sibling_faces = Vec::with_capacity(4);
        let mut current_face = face_dart;
        sibling_faces.push(current_face);

        // Follow the sibling chain to get all 4 faces
        for _ in 1..4 {
            if let Some(sibling_attr) = map.force_read_attribute::<SiblingDartId>(current_face) {
                current_face = sibling_attr.0;
                // Avoid infinite loops if there's a cycle shorter than 4
                if !sibling_faces.contains(&current_face) {
                    sibling_faces.push(current_face);
                } else {
                    break;
                }
            } else {
                break;
            }
        }

        for &sibling in &sibling_faces {
            processed_faces.insert(sibling);
        }

        // Keep the original face_dart as the representative
        sanitized.push(face_dart);
    }

    sanitized
}

fn parallel_refine_with_pairing<T: CoordsFloat>(
    map: &mut CMap2<T>,
    current_level_zipped: impl ParallelIterator<Item = (u32, u32)>,
    next_level: &mut Vec<u32>,
    balance_pile: &mut Vec<u32>,
    geo_verts: &mut [Vertex2<T>],
) {
    let results: Vec<(Vec<u32>, Vec<u32>, Vec<Vertex2<T>>, (usize, usize))> = current_level_zipped
        .map(|(working_dart, start_dart)| {
            let result = refine_with_pairing(map, working_dart, start_dart, geo_verts);

            let valid_children: Vec<u32> = result
                .children
                .into_iter()
                .filter(|&cell| cell != NULL_DART_ID)
                .collect();

            (
                valid_children,
                result.balance_pile,
                result.local_geo_verts,
                result.start_end,
            )
        })
        .collect();

    // Merge results back into shared structures
    for (children, local_pile, local_slice, start_end) in results {
        next_level.extend(children);
        balance_pile.extend(local_pile);
        geo_verts[start_end.0..start_end.1].copy_from_slice(&local_slice);
    }
}

fn balancing<T: CoordsFloat>(
    map: &mut CMap2<T>,
    balance_pile: &mut Vec<u32>,
    geo_verts: &mut [Vertex2<T>],
) {
    while !balance_pile.is_empty() {
        //balance_pile.make_contiguous().shuffle(&mut rand::thread_rng());

        //sanitize the balance pile to avoid processing the same face multiple times including siblings
        let sanitized_balance_pile = sanitize_balance_pile(map, balance_pile);

        let nb_allocations = 64 * sanitized_balance_pile.len();
        let new_dart = map.allocate_unused_darts(nb_allocations);
        let chunks = (0..sanitized_balance_pile.len()).map(|i| new_dart + (i as u32 * 64));

        let balance_pile_zipped: Vec<(u32, u32)> =
            sanitized_balance_pile.into_iter().zip(chunks).collect();

        let balance_pile_zipped = balance_pile_zipped.into_par_iter();

        let results: Vec<(Vec<u32>, Vec<Vertex2<T>>, (usize, usize))> = balance_pile_zipped
            .map(|(face_dart, start_dart)| {
                let res = refine_with_pairing(map, face_dart, start_dart, geo_verts);

                (res.balance_pile, res.local_geo_verts, res.start_end)
            })
            .collect();

        // Merge results back into shared structures
        balance_pile.clear();
        for (local_pile, local_slice, start_end) in results {
            if !local_pile.is_empty() {
                balance_pile.extend(local_pile);
            }
            if start_end != (0, 0) {
                geo_verts[start_end.0..start_end.1].copy_from_slice(&local_slice);
            }
        }
    }
}
