use std::env;

use honeycomb_render::render_2d_map;
use overlay_grid::*;

fn main() {
    let args: Vec<String> = env::args().collect();
    let path = args.get(1).map(|s| s.as_ref()).unwrap_or("manual");
    let map = overlay_grid::<f64>(path, [1., 1.]).unwrap();
    render_2d_map(map);
}
