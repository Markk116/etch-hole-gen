use geo::{Coord, Polygon, Point, Line};
use geo::algorithm::contains::Contains;
use geo::algorithm::euclidean_distance::EuclideanDistance;
use rayon::prelude::*;

/// Generate hexagonal close-packed (HCP) grid points within a polygon
///
/// All coordinates should be in meters (SI base unit).
pub fn generate_hcp_grid(
    polygon: &Polygon<f64>,
    pitch: f64,
    clearance: f64,
) -> Vec<Coord<f64>> {
    // Get bounding box
    let coords: Vec<_> = polygon.exterior().coords().collect();
    if coords.is_empty() {
        return Vec::new();
    }

    let min_x = coords.iter().map(|c| c.x).fold(f64::INFINITY, f64::min);
    let max_x = coords.iter().map(|c| c.x).fold(f64::NEG_INFINITY, f64::max);
    let min_y = coords.iter().map(|c| c.y).fold(f64::INFINITY, f64::min);
    let max_y = coords.iter().map(|c| c.y).fold(f64::NEG_INFINITY, f64::max);

    // HCP grid parameters
    let dx = pitch;
    let dy = pitch * (3.0_f64.sqrt() / 2.0);

    // Generate all candidate positions
    let num_rows = ((max_y - min_y) / dy).ceil() as usize + 1;
    let num_cols = ((max_x - min_x) / dx).ceil() as usize + 1;

    let candidates: Vec<Coord<f64>> = (0..num_rows)
        .flat_map(|row| {
            let y = min_y + (row as f64) * dy;
            let x_offset = if row % 2 == 0 { 0.0 } else { dx / 2.0 };

            (0..num_cols).map(move |col| {
                let x = min_x + x_offset + (col as f64) * dx;
                Coord { x, y }
            })
        })
        .collect();

    // Filter in parallel
    candidates
        .into_par_iter()
        .filter(|coord| {
            if !polygon.contains(coord) {
                return false;
            }
            if clearance > 0.0 {
                distance_to_boundary(coord, polygon) >= clearance
            } else {
                true
            }
        })
        .collect()
}

/// Compute minimum distance from a point to the polygon boundary
fn distance_to_boundary(point: &Coord<f64>, polygon: &Polygon<f64>) -> f64 {
    let pt = Point::from(*point);
    let exterior = polygon.exterior();

    let mut min_dist = f64::INFINITY;

    let coords: Vec<_> = exterior.coords().collect();
    for i in 0..coords.len() - 1 {
        let line = Line::new(*coords[i], *coords[i + 1]);
        let dist = pt.euclidean_distance(&line);
        if dist < min_dist {
            min_dist = dist;
        }
    }

    min_dist
}
