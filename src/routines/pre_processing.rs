//! Step 0 implementation

use std::collections::HashSet;

use honeycomb_core::geometry::CoordsFloat;

use crate::{model::Geometry2, OverlayGridError};

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

//Create a first cell from the overlapping grid size
