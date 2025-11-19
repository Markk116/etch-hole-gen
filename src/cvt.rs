use geo::{Coord, Polygon, LineString, MultiPolygon};
use geo::algorithm::contains::Contains;
use geo::algorithm::intersects::Intersects;
use geo_clipper::Clipper;
use voronator::delaunator::Point as DelaunatorPoint;
use voronator::VoronoiDiagram;
use anyhow::Result;
use crate::svg_output;

/// Statistics from CVT optimization
#[derive(Debug, Clone)]
pub struct CvtStats {
    pub iterations_run: usize,
    pub initial_variance: f64,
    pub final_variance: f64,
    pub final_elongation: f64,
}

/// Compute Centroidal Voronoi Tessellation using Lloyd's algorithm
///
/// All coordinates should be in meters (SI base unit).
pub fn compute_cvt(
    mut points: Vec<Coord<f64>>,
    boundary: &Polygon<f64>,
    max_iterations: usize,
    convergence_threshold: f64,
    debug_svg_prefix: Option<&str>,
) -> Result<(Vec<Coord<f64>>, CvtStats)> {
    let mut variance_history = Vec::new();
    let mut final_elongation = 1.0;
    let mut iterations_run = 0;

    for iter in 0..max_iterations {
        iterations_run = iter + 1;

        // Build Voronoi diagram
        let voronoi_points: Vec<DelaunatorPoint> = points
            .iter()
            .map(|c| DelaunatorPoint { x: c.x, y: c.y })
            .collect();

        let (min_pt, max_pt) = compute_bounding_box(boundary);

        let diagram = match VoronoiDiagram::new(&min_pt, &max_pt, &voronoi_points) {
            Some(d) => d,
            None => break,
        };

        // Compute new centroids
        let mut new_points = Vec::new();
        let mut total_movement = 0.0;
        let mut cell_areas = Vec::new();
        let mut cell_elongations = Vec::new();

        for (i, old_point) in points.iter().enumerate() {
            let cell = get_voronoi_cell(&diagram, i, boundary);

            if let Some(centroid) = compute_polygon_centroid(&cell) {
                if boundary.contains(&centroid) {
                    total_movement += distance(old_point, &centroid);
                    new_points.push(centroid);
                    cell_areas.push(polygon_area(&cell));
                    cell_elongations.push(compute_cell_elongation(&cell, &centroid));
                } else {
                    new_points.push(*old_point);
                    cell_areas.push(polygon_area(&cell));
                    cell_elongations.push(1.0);
                }
            } else {
                new_points.push(*old_point);
                cell_elongations.push(1.0);
            }
        }

        points = new_points;

        // Track metrics
        let variance = compute_variance(&cell_areas);
        variance_history.push(variance);

        final_elongation = if !cell_elongations.is_empty() {
            cell_elongations.iter().sum::<f64>() / cell_elongations.len() as f64
        } else {
            1.0
        };

        // Debug SVG output
        if let Some(prefix) = debug_svg_prefix {
            let debug_path = format!("{}_{:04}.svg", prefix, iter);
            let scale = 1e6; // Scale up for visibility (m to um scale)
            let _ = svg_output::write_voronoi_svg(&debug_path, boundary, &points, scale);
        }

        // Check convergence
        let avg_movement = total_movement / points.len() as f64;
        if avg_movement < convergence_threshold {
            break;
        }
    }

    // Final centroid snap
    let voronoi_points: Vec<DelaunatorPoint> = points
        .iter()
        .map(|c| DelaunatorPoint { x: c.x, y: c.y })
        .collect();

    let (min_pt, max_pt) = compute_bounding_box(boundary);

    if let Some(diagram) = VoronoiDiagram::new(&min_pt, &max_pt, &voronoi_points) {
        let mut final_points = Vec::new();
        let mut elongations = Vec::new();

        for (i, _) in points.iter().enumerate() {
            let cell = get_voronoi_cell(&diagram, i, boundary);

            if let Some(centroid) = compute_polygon_centroid(&cell) {
                if boundary.contains(&centroid) {
                    final_points.push(centroid);
                    elongations.push(compute_cell_elongation(&cell, &centroid));
                } else {
                    final_points.push(points[i]);
                    elongations.push(1.0);
                }
            } else {
                final_points.push(points[i]);
                elongations.push(1.0);
            }
        }

        points = final_points;

        if !elongations.is_empty() {
            final_elongation = elongations.iter().sum::<f64>() / elongations.len() as f64;
        }
    }

    let stats = CvtStats {
        iterations_run,
        initial_variance: *variance_history.first().unwrap_or(&0.0),
        final_variance: *variance_history.last().unwrap_or(&0.0),
        final_elongation,
    };

    Ok((points, stats))
}

fn compute_bounding_box(boundary: &Polygon<f64>) -> (DelaunatorPoint, DelaunatorPoint) {
    let coords: Vec<_> = boundary.exterior().coords().collect();
    let min_pt = DelaunatorPoint {
        x: coords.iter().map(|p| p.x).fold(f64::INFINITY, f64::min),
        y: coords.iter().map(|p| p.y).fold(f64::INFINITY, f64::min),
    };
    let max_pt = DelaunatorPoint {
        x: coords.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max),
        y: coords.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max),
    };
    (min_pt, max_pt)
}

fn get_voronoi_cell(
    diagram: &VoronoiDiagram<DelaunatorPoint>,
    site_index: usize,
    boundary: &Polygon<f64>,
) -> Polygon<f64> {
    let all_cells = diagram.cells();

    if site_index < all_cells.len() {
        let cell = &all_cells[site_index];

        let mut coords: Vec<Coord<f64>> = cell.points()
            .iter()
            .map(|p| Coord { x: p.x, y: p.y })
            .collect();

        if !coords.is_empty() {
            if let Some(&first) = coords.first() {
                coords.push(first);
            }

            let cell_poly = Polygon::new(LineString::new(coords), vec![]);
            return clip_polygon_to_boundary(&cell_poly, boundary);
        }
    }

    Polygon::new(LineString::new(vec![]), vec![])
}

fn clip_polygon_to_boundary(poly: &Polygon<f64>, boundary: &Polygon<f64>) -> Polygon<f64> {
    if !poly.intersects(boundary) {
        return Polygon::new(LineString::new(vec![]), vec![]);
    }

    // Scale factor for geo-clipper (works with integers internally)
    // For meter-scale coordinates (1e-3 to 1e-2), use large scale
    let scale_factor = 1e9;
    let result = poly.intersection(boundary, scale_factor);

    match result {
        MultiPolygon(polys) if !polys.is_empty() => {
            polys.into_iter()
                .max_by(|a, b| {
                    let area_a = polygon_area(a);
                    let area_b = polygon_area(b);
                    area_a.partial_cmp(&area_b).unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap_or_else(|| Polygon::new(LineString::new(vec![]), vec![]))
        }
        _ => poly.clone(),
    }
}

fn compute_polygon_centroid(polygon: &Polygon<f64>) -> Option<Coord<f64>> {
    let coords: Vec<_> = polygon.exterior().coords().collect();

    if coords.len() < 3 {
        return None;
    }

    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut area = 0.0;

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

    if area.abs() < 1e-20 {
        return None;
    }

    cx /= 6.0 * area;
    cy /= 6.0 * area;

    Some(Coord { x: cx, y: cy })
}

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

fn compute_variance(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }

    let mean = values.iter().sum::<f64>() / values.len() as f64;
    values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64
}

fn compute_cell_elongation(cell: &Polygon<f64>, centroid: &Coord<f64>) -> f64 {
    let coords: Vec<_> = cell.exterior().coords().collect();

    if coords.len() < 3 {
        return 1.0;
    }

    let mut min_dist = f64::INFINITY;
    let mut max_dist = 0.0_f64;

    for coord in coords.iter() {
        let dist = distance(centroid, coord);
        if dist > 1e-20 {
            min_dist = min_dist.min(dist);
            max_dist = max_dist.max(dist);
        }
    }

    if min_dist < 1e-20 || min_dist == f64::INFINITY {
        return 1.0;
    }

    max_dist / min_dist
}

fn distance(p1: &Coord<f64>, p2: &Coord<f64>) -> f64 {
    ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
}
