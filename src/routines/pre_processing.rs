//! Step 0 implementation

use std::collections::HashSet;

use honeycomb::core::geometry::CoordsFloat;
use thiserror::Error;
use vtkio::Vtk;

use crate::model::Geometry2;

use crate::model::{GeoVertices, IsIrregular, RefinementLevel, SiblingDartId};
use honeycomb::core::{
    cmap::{CMap2, CMapBuilder},
    geometry::Vertex2,
};

#[derive(Error, Debug)]
/// Errors that can occur during overlay grid operations
pub enum OverlayGridError {
    /// An orientation issue has been detected in the input geometry.
    #[error("boundary isn't consistently oriented - {0}")]
    InconsistentOrientation(&'static str),
    /// The specified geometry does not match one (or more) requirements of the algorithm.
    #[error("input shape isn't conform to requirements - {0}")]
    InvalidShape(&'static str),
    /// The VTK file used to try to build a `Geometry2` object contains invalid data
    /// (per VTK's specification).
    #[error("invalid/corrupted data in the vtk file - {0}")]
    BadVtkData(&'static str),
    /// The VTK file used to try to build a `Geometry2` object contains valid but unsupported data.
    #[error("unsupported data in the vtk file - {0}")]
    UnsupportedVtkData(&'static str),
}

/// Check for orientation issue **per boundary**.
///
/// This function check for the most obvious orientation issue; given a boundary, are all segments making it up
/// oriented consistently. If it is not the case, then there is at least one of:
///
/// - a vertex being the origin of two segment
/// - a vertex being the end-point of two segment
///
/// This does not cover consistent orientation across distinct boundaries (e.g. a geometry with a hole in it).
pub fn detect_orientation_issue<T: CoordsFloat>(
    geometry: &Geometry2<T>,
) -> Result<(), OverlayGridError> {
    let mut origins = HashSet::new();
    let mut endpoints = HashSet::new();

    for (orig, endp) in &geometry.segments {
        if !origins.insert(orig) || !endpoints.insert(endp) {
            return Err(OverlayGridError::InconsistentOrientation(
                "in-boundary inconsistency",
            ));
        }
    }

    Ok(())
}

pub fn compute_overlapping_grid_size<T: CoordsFloat>(
    geometry: &Geometry2<T>,
) -> Result<[T; 4], OverlayGridError> {
    // compute the minimum bounding box
    let (mut min_x, mut max_x, mut min_y, mut max_y): (T, T, T, T) = {
        let Some(tmp) = geometry.vertices.first() else {
            return Err(OverlayGridError::InvalidShape("no vertex in shape"));
        };

        (tmp.x(), tmp.x(), tmp.y(), tmp.y())
    };

    geometry.vertices.iter().for_each(|v| {
        min_x = min_x.min(v.x());
        max_x = max_x.max(v.x()); // may not be optimal
        min_y = min_y.min(v.y()); // don't care
        max_y = max_y.max(v.y());
    });

    if max_x <= min_x {
        return Err(OverlayGridError::InvalidShape(
            "bounding values along X axis are equal",
        ));
    }
    if max_y <= min_y {
        return Err(OverlayGridError::InvalidShape(
            "bounding values along Y axis are equal",
        ));
    }

    Ok([max_x, max_y, min_x, min_y])
}

pub fn manual_grid<T: CoordsFloat>() -> (CMap2<T>, Vec<Vertex2<T>>) {
    let geo_verts = vec![
        Vertex2::<T>(T::from(1.).unwrap(), T::from(1.).unwrap()),
        Vertex2::<T>(T::from(3.).unwrap(), T::from(3.).unwrap()),
        Vertex2::<T>(T::from(1.).unwrap(), T::from(3.).unwrap()),
        Vertex2::<T>(T::from(3.).unwrap(), T::from(1.).unwrap()),
        Vertex2::<T>(T::from(2.).unwrap(), T::from(2.).unwrap()),
    ];

    // --- BUILD THE CMap2
    let map: CMap2<T> = CMapBuilder::<2, _>::from_n_darts(4)
        .add_attribute::<IsIrregular>()
        .add_attribute::<GeoVertices>()
        .add_attribute::<SiblingDartId>()
        .add_attribute::<RefinementLevel>()
        .build()
        .unwrap();

    map.force_link::<1>(1, 2).unwrap();
    map.force_link::<1>(2, 3).unwrap();
    map.force_link::<1>(3, 4).unwrap();
    map.force_link::<1>(4, 1).unwrap();

    map.force_write_vertex(1, Vertex2(T::from(0.0).unwrap(), T::from(0.0).unwrap()));
    map.force_write_vertex(2, Vertex2(T::from(4.0).unwrap(), T::from(0.0).unwrap()));
    map.force_write_vertex(3, Vertex2(T::from(4.0).unwrap(), T::from(4.0).unwrap()));
    map.force_write_vertex(4, Vertex2(T::from(0.0).unwrap(), T::from(4.0).unwrap()));

    map.force_write_attribute::<GeoVertices>(1, GeoVertices((0u32, geo_verts.len() as u32)));
    map.force_write_attribute::<RefinementLevel>(1, RefinementLevel(0));

    (map, geo_verts)
}

pub fn vtk_grid<T: CoordsFloat>(
    file_path: &str,
) -> Result<(CMap2<T>, Vec<Vertex2<T>>), OverlayGridError> {
    // --- IMPORT VTK INPUT
    let geometry_vtk = match Vtk::import(file_path) {
        Ok(vtk) => vtk,
        Err(e) => panic!("E: could not open specified vtk file - {e}"),
    };
    //----/

    // --- BUILD OUR MODEL FROM THE VTK IMPORT
    let geometry = Geometry2::<T>::try_from(geometry_vtk)?;
    //----/

    // --- FIRST DETECTION OF ORIENTATION ISSUES
    detect_orientation_issue(&geometry)?;
    //----/

    // --- Overlapping grid size computation
    let [max_x, max_y, min_x, min_y] = compute_overlapping_grid_size(&geometry)?;
    let _quadtree_bounds = [min_x, min_y, max_x, max_y];
    let width = max_x - min_x;
    let height = max_y - min_y;

    let padding_factor = T::from(0.1).unwrap();
    let x_padding = width * padding_factor;
    let y_padding = height * padding_factor;

    let grid_min_x = min_x - x_padding;
    let grid_max_x = max_x + x_padding;
    let grid_min_y = min_y - y_padding;
    let grid_max_y = max_y + y_padding;

    // Extract geo vertices from the geometry
    let geo_verts = geometry.vertices;

    // --- BUILD THE CMap2
    let map: CMap2<T> = CMapBuilder::<2, _>::from_n_darts(4)
        .add_attribute::<IsIrregular>()
        .add_attribute::<GeoVertices>()
        .add_attribute::<SiblingDartId>()
        .add_attribute::<RefinementLevel>()
        .build()
        .unwrap();

    map.force_link::<1>(1, 2).unwrap();
    map.force_link::<1>(2, 3).unwrap();
    map.force_link::<1>(3, 4).unwrap();
    map.force_link::<1>(4, 1).unwrap();

    // Set vertices based on geometry bounds with padding
    map.force_write_vertex(1, Vertex2(grid_min_x, grid_min_y)); // Bottom-left
    map.force_write_vertex(2, Vertex2(grid_max_x, grid_min_y)); // Bottom-right
    map.force_write_vertex(3, Vertex2(grid_max_x, grid_max_y)); // Top-right
    map.force_write_vertex(4, Vertex2(grid_min_x, grid_max_y)); // Top-left

    map.force_write_attribute::<GeoVertices>(1, GeoVertices((0u32, geo_verts.len() as u32)));
    map.force_write_attribute::<RefinementLevel>(1, RefinementLevel(0));

    Ok((map, geo_verts))
}
