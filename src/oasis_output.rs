use anyhow::Result;
use geo::{Coord, Polygon};
use libreda_db::prelude::*;
use libreda_oasis::OASISStreamWriter;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use voronator::delaunator::Point as DelaunatorPoint;
use voronator::VoronoiDiagram;

/// Meters to nanometers conversion (OASIS typically uses nm)
const M_TO_NM: f64 = 1e9;

/// Maximum segment length in meters (1 micron = 1e-6 meters)
const MAX_SEGMENT_LENGTH: f64 = 1e-6;

/// Calculate number of segments needed for a circle to keep segments under max length
fn calculate_circle_segments(radius: f64) -> usize {
    let circumference = 2.0 * std::f64::consts::PI * radius;
    let min_segments = (circumference / MAX_SEGMENT_LENGTH).ceil() as usize;
    min_segments.max(16) // At least 16 segments for reasonable circle appearance
}

/// Write hole pattern to OASIS file
///
/// Input coordinates should be in meters (SI base unit).
/// Output OASIS will use nanometers as database units.
pub fn write_oasis(
    output_path: &Path,
    boundary: &Polygon<f64>,
    hole_centers: &[Coord<f64>],
    hole_diameter: f64,
    include_outline: bool,
) -> Result<()> {
    // Create a new chip/layout
    let mut chip = Chip::new();

    // Create top cell
    let top_cell = chip.create_cell("TOP".into());

    // Get layer IDs
    let boundary_layer = chip.find_or_create_layer(0, 0);
    let holes_layer = chip.find_or_create_layer(1, 0);

    // Add boundary if requested
    if include_outline {
        let boundary_pts: Vec<Point<i32>> = boundary.exterior().coords()
            .map(|c| Point::new(
                (c.x * M_TO_NM) as i32,
                (c.y * M_TO_NM) as i32
            ))
            .collect();

        if boundary_pts.len() >= 3 {
            let polygon = SimplePolygon::new(boundary_pts);
            chip.insert_shape(&top_cell, &boundary_layer, polygon.into());
        }
    }

    // Add holes as polygons (circles approximated as polygons)
    let hole_radius = hole_diameter / 2.0;
    let num_segments = calculate_circle_segments(hole_radius);

    for center in hole_centers {
        let circle_pts: Vec<Point<i32>> = (0..num_segments)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * (i as f64) / (num_segments as f64);
                let x = center.x + hole_radius * angle.cos();
                let y = center.y + hole_radius * angle.sin();
                Point::new(
                    (x * M_TO_NM) as i32,
                    (y * M_TO_NM) as i32
                )
            })
            .collect();

        let polygon = SimplePolygon::new(circle_pts);
        chip.insert_shape(&top_cell, &holes_layer, polygon.into());
    }

    // Write to file
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    let oasis_writer = OASISStreamWriter::default();
    oasis_writer.write_layout(&mut writer, &chip)
        .map_err(|e| anyhow::anyhow!("Failed to write OASIS: {:?}", e))?;

    Ok(())
}

/// Write Voronoi diagram with point classification to OASIS file
///
/// Input coordinates should be in meters (SI base unit).
/// Output OASIS will use nanometers as database units.
///
/// Layer assignments:
/// - Layer 0: Boundary polygon
/// - Layer 1: Voronoi cell edges
/// - Layer 2: Internal points (blue in SVG)
/// - Layer 3: Rim points (orange in SVG)
/// - Layer 4: Active points (green in SVG)
/// - Layer 5: Isolation region (circle/square)
pub fn write_voronoi_oasis(
    output_path: &Path,
    boundary: &Polygon<f64>,
    hole_centers: &[Coord<f64>],
    isolation_circle: Option<(Coord<f64>, f64)>, // (center, radius)
    isolation_square: Option<(Coord<f64>, f64)>, // (center, side_length)
    internal_points: &[(usize, Coord<f64>)],
    rim_points: &[(usize, Coord<f64>)],
    active_points: &[(usize, Coord<f64>)],
) -> Result<()> {
    // Create a new chip/layout
    let mut chip = Chip::new();

    // Create top cell
    let top_cell = chip.create_cell("TOP".into());

    // Get layer IDs
    let boundary_layer = chip.find_or_create_layer(0, 0);
    let voronoi_layer = chip.find_or_create_layer(1, 0);
    let internal_layer = chip.find_or_create_layer(2, 0);
    let rim_layer = chip.find_or_create_layer(3, 0);
    let active_layer = chip.find_or_create_layer(4, 0);
    let isolation_layer = chip.find_or_create_layer(5, 0);

    // Add boundary polygon
    let boundary_pts: Vec<Point<i32>> = boundary.exterior().coords()
        .map(|c| Point::new(
            (c.x * M_TO_NM) as i32,
            (c.y * M_TO_NM) as i32
        ))
        .collect();

    if boundary_pts.len() >= 3 {
        let polygon = SimplePolygon::new(boundary_pts);
        chip.insert_shape(&top_cell, &boundary_layer, polygon.into());
    }

    // Compute and add Voronoi diagram
    if hole_centers.len() >= 3 {
        let voronoi_points: Vec<DelaunatorPoint> = hole_centers
            .iter()
            .map(|c| DelaunatorPoint { x: c.x, y: c.y })
            .collect();

        let min_pt = DelaunatorPoint {
            x: hole_centers.iter().map(|p| p.x).fold(f64::INFINITY, f64::min),
            y: hole_centers.iter().map(|p| p.y).fold(f64::INFINITY, f64::min),
        };
        let max_pt = DelaunatorPoint {
            x: hole_centers.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max),
            y: hole_centers.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max),
        };

        if let Some(diagram) = VoronoiDiagram::new(&min_pt, &max_pt, &voronoi_points) {
            // Draw Voronoi cells as polygons (like the boundary)
            let all_cells = diagram.cells();
            for cell in all_cells {
                let cell_points = cell.points();
                if cell_points.len() >= 3 {
                    let voronoi_pts: Vec<Point<i32>> = cell_points.iter()
                        .map(|p| Point::new(
                            (p.x * M_TO_NM) as i32,
                            (p.y * M_TO_NM) as i32
                        ))
                        .collect();
                    
                    let polygon = SimplePolygon::new(voronoi_pts);
                    chip.insert_shape(&top_cell, &voronoi_layer, polygon.into());
                }
            }
        }
    }

    // Add isolation region if present
    if let Some((center, radius)) = isolation_circle {
        let num_segments = calculate_circle_segments(radius);
        let circle_pts: Vec<Point<i32>> = (0..num_segments)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * (i as f64) / (num_segments as f64);
                let x = center.x + radius * angle.cos();
                let y = center.y + radius * angle.sin();
                Point::new(
                    (x * M_TO_NM) as i32,
                    (y * M_TO_NM) as i32
                )
            })
            .collect();

        let polygon = SimplePolygon::new(circle_pts);
        chip.insert_shape(&top_cell, &isolation_layer, polygon.into());
    }

    if let Some((center, side_length)) = isolation_square {
        let half_side = side_length / 2.0;
        let square_pts = vec![
            Point::new(
                ((center.x - half_side) * M_TO_NM) as i32,
                ((center.y - half_side) * M_TO_NM) as i32
            ),
            Point::new(
                ((center.x + half_side) * M_TO_NM) as i32,
                ((center.y - half_side) * M_TO_NM) as i32
            ),
            Point::new(
                ((center.x + half_side) * M_TO_NM) as i32,
                ((center.y + half_side) * M_TO_NM) as i32
            ),
            Point::new(
                ((center.x - half_side) * M_TO_NM) as i32,
                ((center.y + half_side) * M_TO_NM) as i32
            ),
        ];

        let polygon = SimplePolygon::new(square_pts);
        chip.insert_shape(&top_cell, &isolation_layer, polygon.into());
    }

    // Helper function to add point markers as small circles
    let add_point_markers = |chip: &mut Chip, 
                             cell, 
                             layer, 
                             points: &[(usize, Coord<f64>)]| {
        // Use a larger radius for point markers so they're visible (1 micrometer)
        let marker_radius = 1.0e-6; // 1 micron in meters
        let num_segments = calculate_circle_segments(marker_radius);
        
        for (_idx, center) in points {
            let circle_pts: Vec<Point<i32>> = (0..num_segments)
                .map(|i| {
                    let angle = 2.0 * std::f64::consts::PI * (i as f64) / (num_segments as f64);
                    let x = center.x + marker_radius * angle.cos();
                    let y = center.y + marker_radius * angle.sin();
                    Point::new(
                        (x * M_TO_NM) as i32,
                        (y * M_TO_NM) as i32
                    )
                })
                .collect();

            let polygon = SimplePolygon::new(circle_pts);
            chip.insert_shape(cell, layer, polygon.into());
        }
    };

    // Add classified points on their respective layers
    let has_classification = !internal_points.is_empty() || !rim_points.is_empty() || !active_points.is_empty();
    
    if has_classification {
        if !internal_points.is_empty() {
            add_point_markers(&mut chip, &top_cell, &internal_layer, internal_points);
        }
        
        if !rim_points.is_empty() {
            add_point_markers(&mut chip, &top_cell, &rim_layer, rim_points);
        }
        
        if !active_points.is_empty() {
            add_point_markers(&mut chip, &top_cell, &active_layer, active_points);
        }
    } else {
        // No classification provided - draw all hole_centers as simple markers
        let marker_radius = 1.0e-6; // 1 micron in meters
        let num_segments = calculate_circle_segments(marker_radius);
        let all_points_layer = chip.find_or_create_layer(2, 0); // Use layer 2 for all points
        
        for center in hole_centers {
            let circle_pts: Vec<Point<i32>> = (0..num_segments)
                .map(|i| {
                    let angle = 2.0 * std::f64::consts::PI * (i as f64) / (num_segments as f64);
                    let x = center.x + marker_radius * angle.cos();
                    let y = center.y + marker_radius * angle.sin();
                    Point::new(
                        (x * M_TO_NM) as i32,
                        (y * M_TO_NM) as i32
                    )
                })
                .collect();

            let polygon = SimplePolygon::new(circle_pts);
            chip.insert_shape(&top_cell, &all_points_layer, polygon.into());
        }
    }

    // Write to file
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    let oasis_writer = OASISStreamWriter::default();
    oasis_writer.write_layout(&mut writer, &chip)
        .map_err(|e| anyhow::anyhow!("Failed to write OASIS: {:?}", e))?;

    Ok(())
}