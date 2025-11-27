use geo::{Coord, Polygon, LineString, MultiPolygon};
use geo::algorithm::contains::Contains;
use geo::algorithm::intersects::Intersects;
use geo_clipper::{Clipper, EndType, JoinType};
use voronator::delaunator::Point as DelaunatorPoint;
use voronator::VoronoiDiagram;
use anyhow::Result;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;
use rayon::prelude::*;
use crate::svg_output;

/// Defines an isolation region where points are frozen during CVT optimization
///
/// Points within the isolation region are classified as either:
/// - Internal points: Points more than 2×pitch from the region boundary (completely excluded from Voronoi computation)
/// - Rim points: Points within 2×pitch of the region boundary (frozen but included in Voronoi to maintain proper tessellation)
///
/// This provides significant performance improvements by reducing the number of points
/// that need to be processed in each CVT iteration.
#[derive(Debug, Clone)]
pub enum IsolationRegion {
    /// Circular isolation region
    Circle { 
        /// Center point of the circle (in meters)
        center: Coord<f64>, 
        /// Radius of the circle (in meters)
        radius: f64 
    },
    /// Square isolation region  
    Square { 
        /// Center point of the square (in meters)
        center: Coord<f64>, 
        /// Side length of the square (in meters)
        side_length: f64 
    }
}

/// Classification of points relative to an isolation region
#[derive(Debug)]
struct PointClassification {
    /// Points fully inside the isolation region (excluded from Voronoi computation)
    internal: Vec<(usize, Coord<f64>)>,
    /// Points near the isolation region boundary (frozen but included in Voronoi)
    rim: Vec<(usize, Coord<f64>)>,
    /// Points outside the isolation region (actively optimized)
    active: Vec<(usize, Coord<f64>)>,
}

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
///
/// # Parameters
///
/// - `points`: Initial point positions to optimize
/// - `boundary`: Bounding polygon that constrains the tessellation
/// - `clearance`: Minimum distance points must maintain from the boundary (in meters). Use 0.0 for no clearance.
/// - `max_iterations`: Maximum number of Lloyd iterations to perform
/// - `convergence_threshold`: Stop when average point movement falls below this value (in meters)
/// - `debug_svg_prefix`: Optional prefix for debug SVG output files
/// - `start_time`: Start time for progress reporting
/// - `interrupted`: Atomic flag to allow early termination
/// - `isolation_region`: Optional region where points are frozen (not moved during optimization)
/// - `pitch`: Characteristic spacing between points (in meters), used to define rim width (2×pitch)
///
/// # Isolation Region Optimization
///
/// When an `isolation_region` is provided, points are classified into three categories:
///
/// 1. Internal points** (> 2×pitch from region boundary): Completely excluded from computation,
///    stored separately and merged back before returning. This provides the main performance benefit.
///
/// 2. Rim points** (within 2×pitch of region boundary): Frozen in place but included in the
///    Voronoi diagram computation to ensure proper tessellation at the region boundary.
///
/// 3. Active points** (outside region): These are the only points that move during optimization.
///
/// The isolation region center is automatically set to the boundary's centroid.
///
/// # Returns
///
/// Returns optimized point positions and statistics about the optimization process.
pub fn compute_cvt(
    points: Vec<Coord<f64>>,
    boundary: &Polygon<f64>,
    clearance: f64,
    max_iterations: usize,
    convergence_threshold: f64,
    debug_svg_prefix: Option<&str>,
    start_time: Instant,
    interrupted: &Arc<AtomicBool>,
    isolation_region: Option<IsolationRegion>,
    pitch: f64,
) -> Result<(Vec<Coord<f64>>, CvtStats)> {
    let mut variance_history = Vec::new();
    let mut final_elongation = 1.0;
    let mut iterations_run = 0;

    // Create shrunk boundary if clearance is specified
    let working_boundary = if clearance > 0.0 {
        // Use Clipper::offset with negative clearance to shrink the boundary inward
        // Scale factor for geo-clipper (works with integers internally)
        let scale_factor = 1e9;
        
        // Clipper::offset parameters:
        // - delta: the offset distance (negative = inward)
        // - join_type: how to join offset lines (Miter for sharp corners)
        // - end_type: how to handle endpoints (ClosedPolygon for polygons)
        // - miter_limit: maximum distance for miter joins (2.0 is reasonable default)
        let offset_result = boundary.offset(
            -clearance,
            JoinType::Miter(1e-9),
            EndType::ClosedPolygon,
            scale_factor,
        );
        
        // The offset operation returns a MultiPolygon, take the largest one
        if let Some(largest) = offset_result.0.into_iter()
            .max_by(|a, b| {
                let area_a = polygon_area(a);
                let area_b = polygon_area(b);
                area_a.partial_cmp(&area_b).unwrap_or(std::cmp::Ordering::Equal)
            }) {
            largest
        } else {
            // If offset produced no polygons, fall back to original boundary
            eprintln!("Warning: clearance offset produced no polygons, using original boundary");
            boundary.clone()
        }
    } else {
        boundary.clone()
    };

    // Classify points if isolation region is provided
    let classification = if let Some(ref region) = isolation_region {
        // Auto-detect center from boundary centroid if needed
        let region_with_center = match region {
            IsolationRegion::Circle { center, radius } => {
                let auto_center = compute_polygon_centroid(&working_boundary).unwrap_or(*center);
                IsolationRegion::Circle { center: auto_center, radius: *radius }
            }
            IsolationRegion::Square { center, side_length } => {
                let auto_center = compute_polygon_centroid(&working_boundary).unwrap_or(*center);
                IsolationRegion::Square { center: auto_center, side_length: *side_length }
            }
        };
        
        Some((classify_points(&points, &region_with_center, pitch), region_with_center))
    } else {
        None
    };

    // Extract working point set and store internal points
    let (mut working_points, internal_points, rim_indices) = if let Some((ref class, _)) = classification {
        // Build working set: rim + active points (preserving their indices for lookup)
        let mut working = Vec::new();
        let mut rim_idx_set = std::collections::HashSet::new();
        
        // Add rim points first
        for (_original_idx, point) in &class.rim {
            working.push(*point);
            rim_idx_set.insert(working.len() - 1);
        }
        
        // Add active points
        for (_original_idx, point) in &class.active {
            working.push(*point);
        }
        
        (working, class.internal.clone(), rim_idx_set)
    } else {
        // No isolation region - all points are active
        (points, Vec::new(), std::collections::HashSet::new())
    };

    // Output initial state before optimization
    if let Some(prefix) = debug_svg_prefix {
        let debug_path = format!("{}_init.svg", prefix);
        let scale = 1e6;
        
        let (iso_circle, iso_square, internal, rim, active) = if let Some((ref class, ref region)) = classification {
            let (circle, square) = match region {
                IsolationRegion::Circle { center, radius } => (Some((*center, *radius)), None),
                IsolationRegion::Square { center, side_length } => (None, Some((*center, *side_length))),
            };
            (circle, square, class.internal.as_slice(), class.rim.as_slice(), class.active.as_slice())
        } else {
            (None, None, &[][..], &[][..], &[][..])
        };
        
        let _ = svg_output::write_voronoi_svg_with_classification(
            &debug_path,
            boundary,
            &working_points,
            scale,
            iso_circle,
            iso_square,
            internal,
            rim,
            active,
        );
    }

    for iter in 0..max_iterations {
        iterations_run = iter + 1;

        // Build Voronoi diagram from working points
        let voronoi_points: Vec<DelaunatorPoint> = working_points
            .iter()
            .map(|c| DelaunatorPoint { x: c.x, y: c.y })
            .collect();

        let (min_pt, max_pt) = compute_bounding_box(&working_boundary);

        let diagram = match VoronoiDiagram::new(&min_pt, &max_pt, &voronoi_points) {
            Some(d) => d,
            None => break,
        };

        // Compute new centroids in parallel, but freeze rim points
        let results: Vec<(Coord<f64>, f64, f64, f64)> = (0..working_points.len())
            .into_par_iter()
            .map(|i| {
                let old_point = working_points[i];
                
                // Check if this is a rim point (frozen)
                if rim_indices.contains(&i) {
                    let cell = get_voronoi_cell(&diagram, i, &working_boundary);
                    let area = polygon_area(&cell);
                    let elongation = if let Some(centroid) = compute_polygon_centroid(&cell) {
                        compute_cell_elongation(&cell, &centroid)
                    } else {
                        1.0
                    };
                    return (old_point, 0.0, area, elongation);
                }
                
                // Active point - compute new position
                let cell = get_voronoi_cell(&diagram, i, &working_boundary);

                if let Some(centroid) = compute_polygon_centroid(&cell) {
                    if working_boundary.contains(&centroid) {
                        let movement = distance(&old_point, &centroid);
                        let area = polygon_area(&cell);
                        let elongation = compute_cell_elongation(&cell, &centroid);
                        (centroid, movement, area, elongation)
                    } else {
                        let area = polygon_area(&cell);
                        (old_point, 0.0, area, 1.0)
                    }
                } else {
                    (old_point, 0.0, 0.0, 1.0)
                }
            })
            .collect();

        // Unpack results
        let mut new_points = Vec::with_capacity(results.len());
        let mut total_movement = 0.0;
        let mut cell_areas = Vec::with_capacity(results.len());
        let mut cell_elongations = Vec::with_capacity(results.len());

        for (point, movement, area, elongation) in results {
            new_points.push(point);
            total_movement += movement;
            if area > 0.0 {
                cell_areas.push(area);
            }
            cell_elongations.push(elongation);
        }

        working_points = new_points;

        // Track metrics
        let variance = compute_variance(&cell_areas);
        variance_history.push(variance);

        final_elongation = if !cell_elongations.is_empty() {
            cell_elongations.iter().sum::<f64>() / cell_elongations.len() as f64
        } else {
            1.0
        };

        // Print progress at 10% intervals
        let print_interval = (max_iterations / 10).max(1);
        if iter % print_interval == 0 || iter == max_iterations - 1 {
            let avg_movement = total_movement / working_points.len() as f64;
            println!("      Iter {}/{}: movement={:.2e}, variance={:.2e} [{:.2}s]",
                iter + 1, max_iterations, avg_movement, variance, start_time.elapsed().as_secs_f64());
        }

        // Debug SVG output
        if let Some(prefix) = debug_svg_prefix {
            let debug_path = format!("{}_{:04}.svg", prefix, iter);
            let scale = 1e6;
            
            let (iso_circle, iso_square, internal, rim, active) = if let Some((ref class, ref region)) = classification {
                let (circle, square) = match region {
                    IsolationRegion::Circle { center, radius } => (Some((*center, *radius)), None),
                    IsolationRegion::Square { center, side_length } => (None, Some((*center, *side_length))),
                };
                (circle, square, class.internal.as_slice(), class.rim.as_slice(), class.active.as_slice())
            } else {
                (None, None, &[][..], &[][..], &[][..])
            };
            
            let _ = svg_output::write_voronoi_svg_with_classification(
                &debug_path,
                boundary,
                &working_points,
                scale,
                iso_circle,
                iso_square,
                internal,
                rim,
                active,
            );
        }

        // Check convergence
        let avg_movement = total_movement / working_points.len() as f64;
        if avg_movement < convergence_threshold {
            break;
        }

        // Check for interrupt
        if interrupted.load(Ordering::SeqCst) {
            println!("      Interrupted after {} iterations", iter + 1);
            break;
        }
    }

    // Merge internal points back with working points
    let final_points = if let Some((ref _class, _)) = classification {
        let mut result = Vec::with_capacity(internal_points.len() + working_points.len());
        
        // Merge all points back (order doesn't matter as long as we don't lose any)
        for (_idx, point) in &internal_points {
            result.push(*point);
        }
        for point in &working_points {
            result.push(*point);
        }
        
        result
    } else {
        working_points
    };

    let stats = CvtStats {
        iterations_run,
        initial_variance: *variance_history.first().unwrap_or(&0.0),
        final_variance: *variance_history.last().unwrap_or(&0.0),
        final_elongation,
    };

    Ok((final_points, stats))
}

pub fn compute_bounding_box(boundary: &Polygon<f64>) -> (DelaunatorPoint, DelaunatorPoint) {
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

/// Check if a point is inside an isolation region
fn point_in_isolation_region(point: &Coord<f64>, region: &IsolationRegion) -> bool {
    match region {
        IsolationRegion::Circle { center, radius } => {
            distance(point, center) <= *radius
        }
        IsolationRegion::Square { center, side_length } => {
            let half_side = side_length / 2.0;
            (point.x - center.x).abs() <= half_side && (point.y - center.y).abs() <= half_side
        }
    }
}

/// Compute distance from a point to the boundary of an isolation region
///
/// Returns positive distance if outside, negative if inside
fn distance_from_isolation_boundary(point: &Coord<f64>, region: &IsolationRegion) -> f64 {
    match region {
        IsolationRegion::Circle { center, radius } => {
            radius - distance(point, center)
        }
        IsolationRegion::Square { center, side_length } => {
            let half_side = side_length / 2.0;
            let dx = (point.x - center.x).abs();
            let dy = (point.y - center.y).abs();
            
            // Distance to nearest edge
            let dist_x = half_side - dx;
            let dist_y = half_side - dy;
            
            dist_x.min(dist_y)
        }
    }
}

/// Classify points into internal, rim, and active sets based on isolation region
///
/// - Internal: Points more than 2×pitch inside the region boundary
/// - Rim: Points within 2×pitch of the region boundary  
/// - Active: Points outside the isolation region
fn classify_points(
    points: &[Coord<f64>],
    region: &IsolationRegion,
    pitch: f64,
) -> PointClassification {
    let rim_width = 2.0 * pitch;
    
    let mut internal = Vec::new();
    let mut rim = Vec::new();
    let mut active = Vec::new();
    
    for (idx, point) in points.iter().enumerate() {
        if point_in_isolation_region(point, region) {
            let dist_from_boundary = distance_from_isolation_boundary(point, region);
            
            if dist_from_boundary > rim_width {
                // Deep inside the isolation region
                internal.push((idx, *point));
            } else {
                // Near the boundary - include in Voronoi but freeze
                rim.push((idx, *point));
            }
        } else {
            // Outside isolation region - actively optimize
            active.push((idx, *point));
        }
    }
    
    PointClassification {
        internal,
        rim,
        active,
    }
}