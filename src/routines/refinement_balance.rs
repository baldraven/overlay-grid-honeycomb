use std::collections::VecDeque;

use crate::{
    model::RefinementLevel,
    routines::{
        edge_splitting::*, helpers::*, pre_processing::*,
    },
};
use honeycomb_core::{
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
    let to_refine = refine_cell(map, 1, geo_verts, &mut balance_pile);
    current_level.push_back(to_refine);

    println!("Starting refinement with 1 initial cell");

    // Process each refinement level
    for depth in 0..max_depth {
        if current_level.is_empty() {
            println!("No more cells to refine at depth {}", depth);
            break;
        };

        println!(
            "Processing refinement level {} with {} cells",
            depth,
            current_level.len()
        );

        let mut next_level = VecDeque::new();

        // Refine all cells at current level
        while let Some(working_dart) = current_level.pop_front() {
            let parent_cells = refine_with_pairing(map, working_dart, geo_verts, &mut balance_pile);

            // Add children to next level
            for child_darts in parent_cells {
                if child_darts != 0 {
                    next_level.push_back(child_darts);
                }
            }
        }

        println!(
            "Level {} complete. {} cells added for next level",
            depth,
            next_level.len()
        );

        balancing(map, &mut balance_pile, geo_verts);

        current_level = next_level;
    }

    if !current_level.is_empty() {
        println!(
            "Refinement stopped at max depth {} with {} cells remaining",
            max_depth,
            current_level.len()
        );
    }

    Ok(())
}

fn balancing<T: CoordsFloat>(
    map: &mut CMap2<T>,
    balance_pile: &mut VecDeque<u32>,
    geo_verts: &mut [Vertex2<T>],
) {
    while !balance_pile.is_empty() {
        let dart = balance_pile.pop_front().unwrap();

        if ensure_not_balanced(map, dart) {
            refine_with_pairing(map, dart, geo_verts, balance_pile);
        } else {
            println!(
                "Skipping balancing for dart {} as it's already scheduled for refinement",
                dart
            );
        }
    }
}

fn ensure_not_balanced<T: CoordsFloat>(map: &CMap2<T>, working_dart: u32) -> bool {
    let face_id = map.face_id(working_dart);
    let face_darts = collect_face_darts(map, working_dart);

    // Get the current cell's refinement level
    let current_level = map
        .force_read_attribute::<RefinementLevel>(face_id)
        .unwrap()
        .0;

    // Check all neighbors to see if any have a refinement level difference > 1
    for &dart in &face_darts {
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

    false // Cell is already balanced
}
