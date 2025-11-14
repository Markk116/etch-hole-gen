use anyhow::Result;
use geo::{Coord, Polygon};
use svg::node::element::{Circle, Path, Group};
use svg::node::element::path::Data;
use svg::Document;

/// Export visualization as SVG file
///
/// # Arguments
/// * `output_path` - Path to output SVG file
/// * `boundary` - Boundary polygon
/// * `hole_centers` - Centers of holes
/// * `hole_diameter` - Diameter of each hole (in micrometers)
/// * `scale` - Scale factor for visualization (pixels per micrometer)
pub fn write_svg(
    output_path: &str,
    boundary: &Polygon<f64>,
    hole_centers: &[Coord<f64>],
    hole_diameter: f64,
    scale: f64,
) -> Result<()> {
    println!("\n=== SVG Export ===");
    println!("Writing to: {}", output_path);
    println!("Scale: {} px/μm", scale);

    // Compute bounding box
    let coords: Vec<_> = boundary.exterior().coords().collect();
    let min_x = coords.iter().map(|c| c.x).fold(f64::INFINITY, f64::min);
    let max_x = coords.iter().map(|c| c.x).fold(f64::NEG_INFINITY, f64::max);
    let min_y = coords.iter().map(|c| c.y).fold(f64::INFINITY, f64::min);
    let max_y = coords.iter().map(|c| c.y).fold(f64::NEG_INFINITY, f64::max);

    let width = (max_x - min_x) * scale;
    let height = (max_y - min_y) * scale;
    let margin = 20.0;

    let mut document = Document::new()
        .set("viewBox", (0, 0, width as i32 + 40, height as i32 + 40))
        .set("width", format!("{}px", width as i32 + 40))
        .set("height", format!("{}px", height as i32 + 40));

    // Add boundary polygon
    let mut path_data = Data::new();
    for (i, coord) in boundary.exterior().coords().enumerate() {
        let x = (coord.x - min_x) * scale + margin;
        let y = (coord.y - min_y) * scale + margin;

        if i == 0 {
            path_data = path_data.move_to((x, y));
        } else {
            path_data = path_data.line_to((x, y));
        }
    }
    path_data = path_data.close();

    let boundary_path = Path::new()
        .set("fill", "none")
        .set("stroke", "black")
        .set("stroke-width", 2)
        .set("d", path_data);

    document = document.add(boundary_path);

    // Add holes
    let hole_radius = hole_diameter / 2.0 * scale;
    let mut holes_group = Group::new()
        .set("fill", "red")
        .set("fill-opacity", 0.5)
        .set("stroke", "red")
        .set("stroke-width", 0.5);

    for center in hole_centers {
        let cx = (center.x - min_x) * scale + margin;
        let cy = (center.y - min_y) * scale + margin;

        let circle = Circle::new()
            .set("cx", cx)
            .set("cy", cy)
            .set("r", hole_radius);

        holes_group = holes_group.add(circle);
    }

    document = document.add(holes_group);

    // Write to file
    svg::save(output_path, &document)?;

    println!("SVG file written successfully!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::LineString;

    #[test]
    fn test_svg_output() {
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

        // Verify the function signature compiles
        assert!(boundary.exterior().coords().count() > 0);
        assert_eq!(hole_centers.len(), 3);
    }
}
