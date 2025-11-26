use anyhow::Result;
use dxf::entities::*;
use dxf::{Drawing, LwPolylineVertex};
use geo::{Coord, Polygon};

/// Meters to millimeters conversion
const M_TO_MM: f64 = 1e3;

/// Write optimized hole pattern to DXF file
///
/// Input coordinates should be in meters (SI base unit).
/// Output DXF will be in millimeters.
pub fn write_dxf(
    output_path: &str,
    boundary: &Polygon<f64>,
    hole_centers: &[Coord<f64>],
    hole_diameter: f64,
    include_outline: bool,
) -> Result<()> {
    let mut drawing = Drawing::new();

    // Optionally add boundary polygon as LWPOLYLINE on layer "BOUNDARY"
    if include_outline {
        add_boundary_to_drawing(&mut drawing, boundary, "BOUNDARY")?;
    }

    // Add holes as circles on layer "HOLES"
    let hole_radius = hole_diameter / 2.0;
    add_holes_to_drawing(&mut drawing, hole_centers, hole_radius, "HOLES")?;

    // Save the drawing
    drawing
        .save_file(output_path)
        .map_err(|e| anyhow::anyhow!("Failed to save DXF: {:?}", e))?;

    Ok(())
}

/// Add boundary polygon to drawing as LWPOLYLINE
fn add_boundary_to_drawing(
    drawing: &mut Drawing,
    boundary: &Polygon<f64>,
    layer: &str,
) -> Result<()> {
    let mut vertices = Vec::new();

    for coord in boundary.exterior().coords() {
        vertices.push(LwPolylineVertex {
            x: coord.x * M_TO_MM,
            y: coord.y * M_TO_MM,
            id: 0,
            starting_width: 0.0,
            ending_width: 0.0,
            bulge: 0.0,
        });
    }

    if !vertices.is_empty() {
        let mut polyline = LwPolyline::default();
        polyline.vertices = vertices;
        polyline.flags = 1; // Closed polyline

        let mut entity = Entity::new(EntityType::LwPolyline(polyline));
        entity.common.layer = layer.to_string();

        drawing.add_entity(entity);
    }

    Ok(())
}

/// Add holes as circles to drawing
fn add_holes_to_drawing(
    drawing: &mut Drawing,
    hole_centers: &[Coord<f64>],
    hole_radius: f64,
    layer: &str,
) -> Result<()> {
    for center in hole_centers {
        let mut circle = Circle::default();
        circle.center = dxf::Point::new(
            center.x * M_TO_MM,
            center.y * M_TO_MM,
            0.0
        );
        circle.radius = hole_radius * M_TO_MM;

        let mut entity = Entity::new(EntityType::Circle(circle));
        entity.common.layer = layer.to_string();

        drawing.add_entity(entity);
    }

    Ok(())
}
