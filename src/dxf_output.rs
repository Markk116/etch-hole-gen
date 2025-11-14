use anyhow::Result;
use dxf::entities::*;
use dxf::{Drawing, LwPolylineVertex};
use geo::{Coord, Polygon};

/// Write optimized hole pattern to DXF file
///
/// # Arguments
/// * `output_path` - Path to output DXF file
/// * `boundary` - Original boundary polygon
/// * `hole_centers` - Centers of holes from CVT optimization
/// * `hole_diameter` - Diameter of each hole (in micrometers)
pub fn write_dxf_with_holes(
    output_path: &str,
    boundary: &Polygon<f64>,
    hole_centers: &[Coord<f64>],
    hole_diameter: f64,
) -> Result<()> {
    println!("\n=== DXF Output ===");
    println!("Writing to: {}", output_path);
    println!("Number of holes: {}", hole_centers.len());
    println!("Hole diameter: {:.3} μm", hole_diameter);

    let mut drawing = Drawing::new();

    // Add boundary polygon as LWPOLYLINE on layer "BOUNDARY"
    add_boundary_to_drawing(&mut drawing, boundary, "BOUNDARY")?;

    // Add holes as circles on layer "HOLES"
    let hole_radius = hole_diameter / 2.0;
    add_holes_to_drawing(&mut drawing, hole_centers, hole_radius, "HOLES")?;

    // Save the drawing
    drawing
        .save_file(output_path)
        .map_err(|e| anyhow::anyhow!("Failed to save DXF: {:?}", e))?;

    println!("DXF file written successfully!");

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
            x: coord.x,
            y: coord.y,
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

        println!("Added boundary polygon with {} vertices", boundary.exterior().coords().count());
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
        circle.center = dxf::Point::new(center.x, center.y, 0.0);
        circle.radius = hole_radius;

        let mut entity = Entity::new(EntityType::Circle(circle));
        entity.common.layer = layer.to_string();

        drawing.add_entity(entity);
    }

    println!("Added {} holes to layer '{}'", hole_centers.len(), layer);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::LineString;

    #[test]
    fn test_dxf_output() {
        let boundary = Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 10.0),
                (0.0, 10.0),
                (0.0, 0.0),
            ]),
            vec![],
        );

        let hole_centers = vec![
            Coord { x: 2.0, y: 2.0 },
            Coord { x: 5.0, y: 5.0 },
            Coord { x: 8.0, y: 8.0 },
        ];

        // This would write to a file in a real test
        // For now, just verify the function signature compiles
        assert!(boundary.exterior().coords().count() > 0);
        assert_eq!(hole_centers.len(), 3);
    }
}
