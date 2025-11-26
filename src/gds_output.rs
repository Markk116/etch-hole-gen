use anyhow::Result;
use geo::{Coord, Polygon};
use gds21::{GdsLibrary, GdsStruct, GdsBoundary, GdsPoint, GdsLayerSpec};
use std::path::Path;

/// Meters to nanometers conversion (GDS typically uses nm as database unit)
const M_TO_NM: f64 = 1e9;

/// Database units per user unit (1 nm)
const DB_UNIT: f64 = 1e-9;
/// User units per meter
const USER_UNIT: f64 = 1e-6;

/// Write hole pattern to GDS file
///
/// Input coordinates should be in meters (SI base unit).
/// Output GDS will use nanometers as database units.
pub fn write_gds(
    output_path: &Path,
    boundary: &Polygon<f64>,
    hole_centers: &[Coord<f64>],
    hole_diameter: f64,
    include_outline: bool,
) -> Result<()> {
    let mut lib = GdsLibrary::new("etch_holes");
    lib.units = gds21::GdsUnits(DB_UNIT, USER_UNIT);

    let mut top_struct = GdsStruct::new("TOP");

    // Layer definitions
    let boundary_layer = GdsLayerSpec { layer: 0, xtype: 0 };
    let holes_layer = GdsLayerSpec { layer: 1, xtype: 0 };

    // Add boundary if requested
    if include_outline {
        let boundary_pts: Vec<GdsPoint> = boundary.exterior().coords()
            .map(|c| GdsPoint::new(
                (c.x * M_TO_NM) as i32,
                (c.y * M_TO_NM) as i32
            ))
            .collect();

        if boundary_pts.len() >= 3 {
            let mut boundary_elem = GdsBoundary::default();
            boundary_elem.layer = boundary_layer.layer;
            boundary_elem.datatype = boundary_layer.xtype;
            boundary_elem.xy = boundary_pts;
            top_struct.elems.push(gds21::GdsElement::GdsBoundary(boundary_elem));
        }
    }

    // Add holes as polygons (circles approximated as polygons)
    let hole_radius = hole_diameter / 2.0;
    let num_segments = 32; // Segments per circle

    for center in hole_centers {
        let circle_pts: Vec<GdsPoint> = (0..=num_segments)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * (i as f64) / (num_segments as f64);
                let x = center.x + hole_radius * angle.cos();
                let y = center.y + hole_radius * angle.sin();
                GdsPoint::new(
                    (x * M_TO_NM) as i32,
                    (y * M_TO_NM) as i32
                )
            })
            .collect();

        let mut hole_elem = GdsBoundary::default();
        hole_elem.layer = holes_layer.layer;
        hole_elem.datatype = holes_layer.xtype;
        hole_elem.xy = circle_pts;
        top_struct.elems.push(gds21::GdsElement::GdsBoundary(hole_elem));
    }

    lib.structs.push(top_struct);

    // Write to file
    lib.save(output_path)
        .map_err(|e| anyhow::anyhow!("Failed to save GDS: {:?}", e))?;

    Ok(())
}
