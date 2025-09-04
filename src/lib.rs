pub(crate) mod model;
pub(crate) mod routines;

use crate::routines::{
        dualization::dualize_map, pre_processing::*,
        refinement_balance::refinement, regularisation::regularize_map,
    };
use honeycomb_core::{
    cmap::CMap2,
    geometry::CoordsFloat,
};

const MAX_DEPTH: u32 = 3;

#[allow(clippy::needless_pass_by_value)]
/// Create an overlay grid from a VTK file input
pub fn overlay_grid<T: CoordsFloat>(
    file_path: impl AsRef<std::path::Path>,
    _grid_cell_sizes: [T; 2],
) -> Result<CMap2<T>, OverlayGridError> {
    // Determine whether to use VTK file or manual grid based on file path
    let file_path_str = file_path.as_ref().to_string_lossy();
    let (mut map, mut geo_verts) = if file_path_str == "manual" || file_path_str.is_empty() {
        manual_grid::<T>()
    } else {
        vtk_grid::<T>(&file_path_str)?
    };

    refinement(&mut map, &mut geo_verts, MAX_DEPTH)?;
    regularize_map(&mut map);
    let dual_map = dualize_map(&map);

    Ok(dual_map)
}

#[cfg(test)]
mod test;
