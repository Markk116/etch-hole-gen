use geo::{Coord, Polygon, LineString, MultiPolygon};
use geo::algorithm::contains::Contains;
use geo::algorithm::intersects::Intersects;
use geo_clipper::Clipper;
use voronator::delaunator::Point as DelaunatorPoint;
use voronator::VoronoiDiagram;
use anyhow::Result;
use crate::svg_output;

/// Compute Centroidal Voronoi Tessellation using Lloyd's algorithm
///
/// # Arguments
/// * `initial_points` - Starting point distribution (e.g., from HCP grid)
/// * `boundary` - Polygon boundary to constrain the tessellation
/// * `max_iterations` - Maximum number of Lloyd iterations
/// * `convergence_threshold` - Stop when centroid movement is below this threshold
/// * `debug_svg_prefix` - Optional prefix for debug SVG output at each iteration
///
/// # Returns
/// Optimized point positions and convergence metrics
pub fn compute_cvt(
    mut points: Vec<Coord<f64>>,
    boundary: &Polygon<f64>,
    max_iterations: usize,
    convergence_threshold: f64,
    debug_svg_prefix: Option<&str>,
) -> Result<(Vec<Coord<f64>>, Vec<f64>)> {
    println!("\n=== CVT Computation (Lloyd's Algorithm) ===");
    println!("Initial points: {}", points.len());
    println!("Max iterations: {}", max_iterations);
    println!("Convergence threshold: {:.6}", convergence_threshold);

    let mut variance_history = Vec::new();
    let mut elongation_history = Vec::new();

    for iter in 0..max_iterations {
        // Convert points to voronator format
        let voronoi_points: Vec<DelaunatorPoint> = points
            .iter()
            .map(|c| DelaunatorPoint { x: c.x, y: c.y })
            .collect();

        // Compute bounding box for Voronoi diagram
        let min_pt = DelaunatorPoint {
            x: points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min),
            y: points.iter().map(|p| p.y).fold(f64::INFINITY, f64::min),
        };
        let max_pt = DelaunatorPoint {
            x: points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max),
            y: points.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max),
        };

        // Compute Voronoi diagram
        let diagram = match VoronoiDiagram::new(&min_pt, &max_pt, &voronoi_points) {
            Some(d) => d,
            None => {
                eprintln!("Voronoi computation failed at iteration {}", iter);
                break;
            }
        };

        // For each point, compute the centroid of its Voronoi cell clipped to boundary
        let mut new_points = Vec::new();
        let mut total_movement = 0.0;
        let mut cell_areas = Vec::new();
        let mut cell_elongations = Vec::new();

        for (i, old_point) in points.iter().enumerate() {
            // Get Voronoi cell for this point
            let cell = get_voronoi_cell(&diagram, i, boundary);

            // Compute centroid of the clipped cell
            if let Some(centroid) = compute_polygon_centroid(&cell) {
                // Ensure centroid is within boundary
                if boundary.contains(&centroid) {
                    let movement = distance(old_point, &centroid);
                    total_movement += movement;
                    new_points.push(centroid);

                    // Track cell area for variance computation
                    let area = polygon_area(&cell);
                    cell_areas.push(area);

                    // Track cell elongation (measure of circularity)
                    let elongation = compute_cell_elongation(&cell, &centroid);
                    cell_elongations.push(elongation);
                } else {
                    // Keep original point if centroid is outside
                    new_points.push(*old_point);
                    let area = polygon_area(&cell);
                    cell_areas.push(area);
                    cell_elongations.push(1.0);
                }
            } else {
                // Keep original point if centroid computation fails
                new_points.push(*old_point);
                cell_elongations.push(1.0);
            }
        }

        points = new_points;

        // Compute variance of cell areas
        let variance = compute_variance(&cell_areas);
        variance_history.push(variance);

        // Compute average elongation (1.0 = circular, higher = more elongated)
        let avg_elongation = if !cell_elongations.is_empty() {
            cell_elongations.iter().sum::<f64>() / cell_elongations.len() as f64
        } else {
            1.0
        };
        elongation_history.push(avg_elongation);

        let avg_movement = total_movement / points.len() as f64;

        if iter % 5 == 0 || iter == max_iterations - 1 {
            println!(
                "Iteration {}: avg movement = {:.6}, variance = {:.6}, elongation = {:.3}",
                iter, avg_movement, variance, avg_elongation
            );
        }

        // Output debug SVG if requested
        if let Some(prefix) = debug_svg_prefix {
            let debug_path = format!("{}_{:04}.svg", prefix, iter);
            if let Err(e) = svg_output::write_voronoi_svg(&debug_path, boundary, &points, 10.0) {
                eprintln!("Warning: Failed to write debug SVG {}: {}", debug_path, e);
            }
        }

        // Check convergence
        if avg_movement < convergence_threshold {
            println!(
                "\nConverged after {} iterations (avg movement: {:.6})",
                iter, avg_movement
            );
            break;
        }
    }

    // Final step: Snap points to exact centroids of their Voronoi cells
    println!("\n=== Final Centroid Positioning ===");
    let voronoi_points: Vec<DelaunatorPoint> = points
        .iter()
        .map(|c| DelaunatorPoint { x: c.x, y: c.y })
        .collect();

    let min_pt = DelaunatorPoint {
        x: points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min),
        y: points.iter().map(|p| p.y).fold(f64::INFINITY, f64::min),
    };
    let max_pt = DelaunatorPoint {
        x: points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max),
        y: points.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max),
    };

    if let Some(diagram) = VoronoiDiagram::new(&min_pt, &max_pt, &voronoi_points) {
        let mut final_points = Vec::new();
        let mut final_elongations = Vec::new();

        for (i, _) in points.iter().enumerate() {
            let cell = get_voronoi_cell(&diagram, i, boundary);

            if let Some(centroid) = compute_polygon_centroid(&cell) {
                if boundary.contains(&centroid) {
                    final_points.push(centroid);
                    let elongation = compute_cell_elongation(&cell, &centroid);
                    final_elongations.push(elongation);
                } else {
                    final_points.push(points[i]);
                    final_elongations.push(1.0);
                }
            } else {
                final_points.push(points[i]);
                final_elongations.push(1.0);
            }
        }

        points = final_points;

        let avg_final_elongation = if !final_elongations.is_empty() {
            final_elongations.iter().sum::<f64>() / final_elongations.len() as f64
        } else {
            1.0
        };

        println!("Final average elongation: {:.3}", avg_final_elongation);
    }

    println!("Final points: {}", points.len());

    Ok((points, variance_history))
}

/// Extract Voronoi cell for a specific site and clip to boundary
fn get_voronoi_cell(
    diagram: &VoronoiDiagram<DelaunatorPoint>,
    site_index: usize,
    boundary: &Polygon<f64>,
) -> Polygon<f64> {
    // Get all cells and access the one we need
    let all_cells = diagram.cells();

    if site_index < all_cells.len() {
        let cell = &all_cells[site_index];

        // Convert voronator polygon to geo polygon
        let mut coords: Vec<Coord<f64>> = cell.points()
            .iter()
            .map(|p| Coord { x: p.x, y: p.y })
            .collect();

        if !coords.is_empty() {
            // Close the polygon
            if let Some(&first) = coords.first() {
                coords.push(first);
            }

            let cell_poly = Polygon::new(LineString::new(coords), vec![]);

            // Clip the Voronoi cell to the boundary polygon
            return clip_polygon_to_boundary(&cell_poly, boundary);
        }
    }

    // Return empty polygon if cell not found
    Polygon::new(LineString::new(vec![]), vec![])
}

/// Clip a polygon to the boundary using intersection
fn clip_polygon_to_boundary(poly: &Polygon<f64>, boundary: &Polygon<f64>) -> Polygon<f64> {
    // Quick check: if poly doesn't intersect boundary at all, return empty
    if !poly.intersects(boundary) {
        return Polygon::new(LineString::new(vec![]), vec![]);
    }

    // Use geo-clipper to compute intersection
    // Scale factor: clipper works with integer coordinates, so we scale up
    // For mm-scale coordinates (0-4mm), use larger scale factor
    let scale_factor = 1000.0; // 1000x scale for better precision
    let result = poly.intersection(boundary, scale_factor);

    // Extract the first polygon from the result (there should typically be only one)
    match result {
        MultiPolygon(polys) if !polys.is_empty() => {
            // Return the largest polygon if there are multiple
            polys.into_iter()
                .max_by(|a, b| {
                    let area_a = polygon_area(a);
                    let area_b = polygon_area(b);
                    area_a.partial_cmp(&area_b).unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap_or_else(|| Polygon::new(LineString::new(vec![]), vec![]))
        }
        _ => {
            // Clipping failed - return unclipped cell as fallback
            poly.clone()
        }
    }
}

/// Compute centroid of a polygon
fn compute_polygon_centroid(polygon: &Polygon<f64>) -> Option<Coord<f64>> {
    let coords: Vec<_> = polygon.exterior().coords().collect();

    if coords.len() < 3 {
        return None;
    }

    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut area = 0.0;

    // Use the formula for polygon centroid
    for i in 0..coords.len() - 1 {
        let x0 = coords[i].x;
        let y0 = coords[i].y;
        let x1 = coords[i + 1].x;
        let y1 = coords[i + 1].y;

        let cross = x0 * y1 - x1 * y0;
        area += cross;
        cx += (x0 + x1) * cross;
        cy += (y0 + y1) * cross;
    }

    area *= 0.5;

    if area.abs() < 1e-10 {
        return None;
    }

    cx /= 6.0 * area;
    cy /= 6.0 * area;

    Some(Coord { x: cx, y: cy })
}

/// Compute area of a polygon
fn polygon_area(polygon: &Polygon<f64>) -> f64 {
    let coords: Vec<_> = polygon.exterior().coords().collect();

    if coords.len() < 3 {
        return 0.0;
    }

    let mut area = 0.0;

    for i in 0..coords.len() - 1 {
        area += coords[i].x * coords[i + 1].y - coords[i + 1].x * coords[i].y;
    }

    (area * 0.5).abs()
}

/// Compute variance of a set of values
fn compute_variance(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }

    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64;

    variance
}

/// Compute cell elongation (measure of circularity)
///
/// Elongation is the ratio of max distance to min distance from centroid to cell vertices.
/// - 1.0 = perfect circle (all vertices equidistant)
/// - >1.0 = elongated (higher = more elongated)
fn compute_cell_elongation(cell: &Polygon<f64>, centroid: &Coord<f64>) -> f64 {
    let coords: Vec<_> = cell.exterior().coords().collect();

    if coords.len() < 3 {
        return 1.0;
    }

    let mut min_dist = f64::INFINITY;
    let mut max_dist = 0.0_f64;

    for coord in coords.iter() {
        let dist = distance(centroid, coord);
        if dist > 1e-10 {  // Avoid zero distances
            min_dist = min_dist.min(dist);
            max_dist = max_dist.max(dist);
        }
    }

    if min_dist < 1e-10 || min_dist == f64::INFINITY {
        return 1.0;
    }

    max_dist / min_dist
}

/// Euclidean distance between two points
fn distance(p1: &Coord<f64>, p2: &Coord<f64>) -> f64 {
    ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polygon_centroid() {
        // Unit square centered at origin
        let square = Polygon::new(
            LineString::from(vec![
                (-1.0, -1.0),
                (1.0, -1.0),
                (1.0, 1.0),
                (-1.0, 1.0),
                (-1.0, -1.0),
            ]),
            vec![],
        );

        let centroid = compute_polygon_centroid(&square).unwrap();
        assert!((centroid.x).abs() < 1e-10);
        assert!((centroid.y).abs() < 1e-10);
    }

    #[test]
    fn test_polygon_area() {
        // Unit square
        let square = Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (1.0, 0.0),
                (1.0, 1.0),
                (0.0, 1.0),
                (0.0, 0.0),
            ]),
            vec![],
        );

        let area = polygon_area(&square);
        assert!((area - 1.0).abs() < 1e-10);
    }
}
