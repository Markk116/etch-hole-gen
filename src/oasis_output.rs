use anyhow::Result;
use geo::{Coord, Polygon};
use libreda_db::prelude::*;
use libreda_oasis::OASISStreamWriter;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

/// Meters to nanometers conversion (OASIS typically uses nm)
const M_TO_NM: f64 = 1e9;

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
    let num_segments = 32;

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
