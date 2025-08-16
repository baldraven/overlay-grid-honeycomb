use std::env;

use honeycomb_render::render_2d_map;
use overlay_grid::*;

fn main() {
    let args: Vec<String> = env::args().collect();

    if let Some(path) = args.get(1) {
        let map = overlay_grid::<f64>(path, [1., 1.]).unwrap();

        render_2d_map(map);
    } else {
        println!(
            "No input geometry specified - you can pass a path to a vtk input as command line argument"
        )
    }
}
