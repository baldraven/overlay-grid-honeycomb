//! Data model for overlay grid operations
//!
//! This module contains data structures for representing geometry and related operations.

use honeycomb_core::{
    cmap::{DartIdType, FaceIdType},
    geometry::{CoordsFloat, Vertex2},
};
use vtkio::{
    model::{CellType, DataSet, VertexNumbers},
    IOBuffer, Vtk,
};

use crate::OverlayGridError;

/// Geometry representation structure.
///
/// For specification of the accepted VTK file format, see [`crate`]'s documentation entry.
pub struct Geometry2<T: CoordsFloat> {
    /// Vertices of the geometry.
    pub vertices: Vec<Vertex2<T>>,
    /// Edges / segments making up the geometry.
    pub segments: Vec<(usize, usize)>,
}

macro_rules! build_vertices {
    ($v: ident) => {{
        if $v.len() % 3 != 0 {
            return Err(OverlayGridError::BadVtkData(
                "vertex list contains an incomplete tuple",
            ));
        }
        $v.chunks_exact(3)
            .map(|slice| {
                // WE IGNORE Z values
                let &[x, y, _] = slice else { unreachable!() };
                Vertex2::from((T::from(x).unwrap(), T::from(y).unwrap()))
            })
            .collect()
    }};
}

/// For specification of the accepted VTK file format, see [`crate`]'s documentation entry.
impl<T: CoordsFloat> TryFrom<Vtk> for Geometry2<T> {
    type Error = OverlayGridError;

    #[allow(clippy::too_many_lines)]
    fn try_from(value: Vtk) -> Result<Self, Self::Error> {
        // What we are reading / how we construct the geometry:
        // The input VTK file should describe boundaries (e.g. edges in 2D) & key vertices (e.g. sharp corners)
        // Those should be described by using simple
        match value.data {
            DataSet::ImageData { .. }
            | DataSet::StructuredGrid { .. }
            | DataSet::RectilinearGrid { .. }
            | DataSet::Field { .. }
            | DataSet::PolyData { .. } => Err(OverlayGridError::UnsupportedVtkData(
                "dataset not supported",
            )),
            DataSet::UnstructuredGrid { pieces, .. } => {
                let mut vertices = Vec::new();
                let mut segments = Vec::new();
                let tmp = pieces.iter().map(|piece| {
                    // assume inline data
                    let Ok(tmp) = piece.load_piece_data(None) else {
                        return Err(OverlayGridError::UnsupportedVtkData(
                            "not inlined data piece",
                        ));
                    };

                    // build vertex list
                    // since we're expecting coordinates, we'll assume floating type
                    // we're also converting directly to our vertex type since we're building a 2-map
                    let vertices: Vec<Vertex2<T>> = match tmp.points {
                        IOBuffer::F64(v) => build_vertices!(v),
                        IOBuffer::F32(v) => build_vertices!(v),
                        _ => {
                            return Err(OverlayGridError::UnsupportedVtkData(
                                "not float or double coordinate representation type",
                            ));
                        }
                    };
                    let mut poi: Vec<usize> = Vec::new();
                    let mut segments: Vec<(usize, usize)> = Vec::new();

                    let vtkio::model::Cells { cell_verts, types } = tmp.cells;
                    match cell_verts {
                        VertexNumbers::Legacy {
                            num_cells,
                            vertices: verts,
                        } => {
                            // check basic stuff
                            if num_cells as usize != types.len() {
                                return Err(OverlayGridError::BadVtkData(
                                    "different # of cells in CELLS and CELL_TYPES",
                                ));
                            }

                            // build a collection of vertex lists corresponding of each cell
                            let mut cell_components: Vec<Vec<usize>> = Vec::new();
                            let mut take_next = 0;
                            for vertex_id in &verts {
                                if take_next == 0 {
                                    // making it usize since it's a counter
                                    take_next = *vertex_id as usize;
                                    cell_components.push(Vec::with_capacity(take_next));
                                } else {
                                    cell_components
                                        .last_mut()
                                        .expect("E: unreachable")
                                        .push(*vertex_id as usize);
                                    take_next -= 1;
                                }
                            }
                            assert_eq!(num_cells as usize, cell_components.len());

                            if let Some(err) = types.iter().zip(cell_components.iter()).find_map(
                                |(cell_type, vids)| match cell_type {
                                    CellType::Vertex => {
                                        if vids.len() != 1 {
                                            return Some(OverlayGridError::BadVtkData(
                                                "`Vertex` with incorrect # of vertices (!=1)",
                                            ));
                                        }
                                        poi.push(vids[0]);
                                        None
                                    }
                                    CellType::PolyVertex => {
                                        Some(OverlayGridError::UnsupportedVtkData(
                                            "`PolyVertex` cell type",
                                        ))
                                    }
                                    CellType::Line => {
                                        if vids.len() != 2 {
                                            return Some(OverlayGridError::BadVtkData(
                                                "`Line` with incorrect # of vertices (!=2)",
                                            ));
                                        }
                                        segments.push((vids[0], vids[1]));
                                        None
                                    }
                                    CellType::PolyLine => {
                                        Some(OverlayGridError::BadVtkData("`PolyLine` cell type"))
                                    }
                                    _ => None, // silent ignore all other cells that do not make up boundaries
                                },
                            ) {
                                return Err(err);
                            }
                        }
                        VertexNumbers::XML { .. } => {
                            return Err(OverlayGridError::UnsupportedVtkData("XML format"));
                        }
                    }
                    Ok((vertices, segments))
                });

                if let Some(e) = tmp.clone().find(Result::is_err) {
                    return Err(e.unwrap_err());
                }

                tmp.filter_map(Result::ok).for_each(|(mut ver, mut seg)| {
                    vertices.append(&mut ver);
                    segments.append(&mut seg);
                });

                Ok(Geometry2 { vertices, segments })
            }
        }
    }
}

use honeycomb_core::attributes::{AttrSparseVec, AttributeBind, AttributeError, AttributeUpdate};
use honeycomb_core::cmap::OrbitPolicy;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct IsIrregular(pub bool);

impl AttributeUpdate for IsIrregular {
    fn merge(attr1: Self, attr2: Self) -> Result<Self, AttributeError> {
        Ok(attr1)
    }

    fn split(attr: Self) -> Result<(Self, Self), AttributeError> {
        Ok((attr, attr))
    }
}

impl AttributeBind for IsIrregular {
    type StorageType = AttrSparseVec<Self>;
    type IdentifierType = DartIdType;
    const BIND_POLICY: OrbitPolicy = OrbitPolicy::Edge;
}

#[derive(Debug, Clone, Default, Copy, PartialEq)]
pub struct GeoVertices(pub (u32, u32));

impl AttributeUpdate for GeoVertices {
    fn merge(attr1: Self, attr2: Self) -> Result<Self, AttributeError> {
        Ok(attr1)
    }

    fn split(attr: Self) -> Result<(Self, Self), AttributeError> {
        Ok((attr, attr))
    }
}

impl AttributeBind for GeoVertices {
    type StorageType = AttrSparseVec<Self>;
    type IdentifierType = FaceIdType;
    const BIND_POLICY: OrbitPolicy = OrbitPolicy::Face;
}

#[derive(Debug, Clone, Default, Copy, PartialEq)]
pub struct SiblingDartId(pub u32);

impl AttributeUpdate for SiblingDartId {
    fn merge(attr1: Self, attr2: Self) -> Result<Self, AttributeError> {
        Ok(attr1)
    }

    fn split(attr: Self) -> Result<(Self, Self), AttributeError> {
        Ok((attr, attr))
    }
}

impl AttributeBind for SiblingDartId {
    type StorageType = AttrSparseVec<Self>;
    type IdentifierType = FaceIdType;
    const BIND_POLICY: OrbitPolicy = OrbitPolicy::Face;
}

#[derive(Debug, Clone, Default, Copy, PartialEq)]
pub struct RefinementLevel(pub u16);

impl AttributeUpdate for RefinementLevel {
    fn merge(attr1: Self, attr2: Self) -> Result<Self, AttributeError> {
        Ok(attr1)
    }

    fn split(attr: Self) -> Result<(Self, Self), AttributeError> {
        Ok((attr, attr))
    }
}

impl AttributeBind for RefinementLevel {
    type StorageType = AttrSparseVec<Self>;
    type IdentifierType = FaceIdType;
    const BIND_POLICY: OrbitPolicy = OrbitPolicy::Face;
}
