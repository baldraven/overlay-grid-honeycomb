use clap::Parser;
use honeycomb::render::render_2d_map;
use overlay_grid::*;

#[derive(Parser)]
#[command(name = "overlay-grid-honeycomb")]
#[command(about = "Overlay grid honeycomb application")]
struct Args {
    /// VTK file path to use as input
    #[arg(long, short = 'f', conflicts_with = "nb_verts")]
    vtk_file: Option<String>,
    
    /// Number of vertices to generate for random grid
    #[arg(long, short = 'n', conflicts_with = "vtk_file", default_value = "100")]
    nb_verts: Option<usize>,
    
    /// Maximum refinement depth
    #[arg(long, short = 'd', default_value = "10")]
    max_depth: u32,
}

fn main() {
    bind_threads();
    let args = Args::parse();
    
    let (path, nb_verts) = if let Some(vtk_path) = args.vtk_file {
        (vtk_path, None)
    } else {
        ("random".to_string(), args.nb_verts)
    };
    
    let map = overlay_grid::<f64>(path, [1., 1.], nb_verts, Some(args.max_depth)).unwrap();
    render_2d_map(map);
}
