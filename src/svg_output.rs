use anyhow::Result;
use geo::{Coord, Polygon};
use svg::node::element::{Circle, Path, Group, Rectangle};
use svg::node::element::path::Data;
use svg::Document;
use voronator::delaunator::Point as DelaunatorPoint;
use voronator::VoronoiDiagram;

/// Export visualization as SVG file
///
/// Input coordinates should be in meters (SI base unit).
/// Scale converts meters to pixels.
pub fn write_svg(
    output_path: &str,
    boundary: &Polygon<f64>,
    hole_centers: &[Coord<f64>],
    hole_diameter: f64,
    scale: f64,
) -> Result<()> {
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

    svg::save(output_path, &document)?;

    Ok(())
}

/// Export Voronoi diagram visualization as SVG file
///
/// Input coordinates should be in meters (SI base unit).
/// Scale converts meters to pixels.
pub fn write_voronoi_svg(
    output_path: &str,
    boundary: &Polygon<f64>,
    hole_centers: &[Coord<f64>],
    scale: f64,
) -> Result<()> {
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

    // Compute Voronoi diagram
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
        // Draw Voronoi cells
        let mut cells_group = Group::new()
            .set("fill", "none")
            .set("stroke", "blue")
            .set("stroke-width", 0.5)
            .set("stroke-opacity", 0.5);

        let all_cells = diagram.cells();
        for cell in all_cells {
            let cell_points = cell.points();
            if cell_points.len() >= 3 {
                let mut path_data = Data::new();
                for (i, p) in cell_points.iter().enumerate() {
                    let x = (p.x - min_x) * scale + margin;
                    let y = (p.y - min_y) * scale + margin;

                    if i == 0 {
                        path_data = path_data.move_to((x, y));
                    } else {
                        path_data = path_data.line_to((x, y));
                    }
                }
                path_data = path_data.close();

                let cell_path = Path::new().set("d", path_data);
                cells_group = cells_group.add(cell_path);
            }
        }
        document = document.add(cells_group);
    }

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

    // Add site points
    let mut sites_group = Group::new()
        .set("fill", "red")
        .set("stroke", "none");

    for center in hole_centers {
        let cx = (center.x - min_x) * scale + margin;
        let cy = (center.y - min_y) * scale + margin;

        let circle = Circle::new()
            .set("cx", cx)
            .set("cy", cy)
            .set("r", 2.0);

        sites_group = sites_group.add(circle);
    }

    document = document.add(sites_group);

    svg::save(output_path, &document)?;

    Ok(())
}

/// Export Voronoi diagram visualization with optional isolation region
///
/// This enhanced version can visualize:
/// - Voronoi cells
/// - Point classification (internal/blue, rim/orange, active/green)
/// - Isolation region boundary (dashed circle or square)
///
/// Input coordinates should be in meters (SI base unit).
/// Scale converts meters to pixels.
pub fn write_voronoi_svg_with_classification(
    output_path: &str,
    boundary: &Polygon<f64>,
    hole_centers: &[Coord<f64>],
    scale: f64,
    isolation_circle: Option<(Coord<f64>, f64)>, // (center, radius)
    isolation_square: Option<(Coord<f64>, f64)>, // (center, side_length)
    internal_points: &[(usize, Coord<f64>)],
    rim_points: &[(usize, Coord<f64>)],
    active_points: &[(usize, Coord<f64>)],
) -> Result<()> {
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

    // Compute Voronoi diagram
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
        // Draw Voronoi cells
        let mut cells_group = Group::new()
            .set("fill", "none")
            .set("stroke", "blue")
            .set("stroke-width", 0.5)
            .set("stroke-opacity", 0.5);

        let all_cells = diagram.cells();
        for cell in all_cells {
            let cell_points = cell.points();
            if cell_points.len() >= 3 {
                let mut path_data = Data::new();
                for (i, p) in cell_points.iter().enumerate() {
                    let x = (p.x - min_x) * scale + margin;
                    let y = (p.y - min_y) * scale + margin;

                    if i == 0 {
                        path_data = path_data.move_to((x, y));
                    } else {
                        path_data = path_data.line_to((x, y));
                    }
                }
                path_data = path_data.close();

                let cell_path = Path::new().set("d", path_data);
                cells_group = cells_group.add(cell_path);
            }
        }
        document = document.add(cells_group);
    }

    // Draw isolation region if present
    if let Some((center, radius)) = isolation_circle {
        let cx = (center.x - min_x) * scale + margin;
        let cy = (center.y - min_y) * scale + margin;
        let r = radius * scale;

        let circle = Circle::new()
            .set("cx", cx)
            .set("cy", cy)
            .set("r", r)
            .set("fill", "none")
            .set("stroke", "purple")
            .set("stroke-width", 2)
            .set("stroke-dasharray", "5,5")
            .set("stroke-opacity", 0.7);

        document = document.add(circle);
    }

    if let Some((center, side_length)) = isolation_square {
        let half_side = side_length / 2.0;
        let x = (center.x - half_side - min_x) * scale + margin;
        let y = (center.y - half_side - min_y) * scale + margin;
        let size = side_length * scale;

        let rect = Rectangle::new()
            .set("x", x)
            .set("y", y)
            .set("width", size)
            .set("height", size)
            .set("fill", "none")
            .set("stroke", "purple")
            .set("stroke-width", 2)
            .set("stroke-dasharray", "5,5")
            .set("stroke-opacity", 0.7);

        document = document.add(rect);
    }

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

    // Add site points with color coding based on classification
    if !internal_points.is_empty() || !rim_points.is_empty() || !active_points.is_empty() {
        // Internal points (blue)
        if !internal_points.is_empty() {
            let mut internal_group = Group::new()
                .set("fill", "blue")
                .set("stroke", "none");

            for (_idx, center) in internal_points {
                let cx = (center.x - min_x) * scale + margin;
                let cy = (center.y - min_y) * scale + margin;

                let circle = Circle::new()
                    .set("cx", cx)
                    .set("cy", cy)
                    .set("r", 2.5);

                internal_group = internal_group.add(circle);
            }
            document = document.add(internal_group);
        }

        // Rim points (orange)
        if !rim_points.is_empty() {
            let mut rim_group = Group::new()
                .set("fill", "orange")
                .set("stroke", "none");

            for (_idx, center) in rim_points {
                let cx = (center.x - min_x) * scale + margin;
                let cy = (center.y - min_y) * scale + margin;

                let circle = Circle::new()
                    .set("cx", cx)
                    .set("cy", cy)
                    .set("r", 2.5);

                rim_group = rim_group.add(circle);
            }
            document = document.add(rim_group);
        }

        // Active points (green)
        if !active_points.is_empty() {
            let mut active_group = Group::new()
                .set("fill", "green")
                .set("stroke", "none");

            for (_idx, center) in active_points {
                let cx = (center.x - min_x) * scale + margin;
                let cy = (center.y - min_y) * scale + margin;

                let circle = Circle::new()
                    .set("cx", cx)
                    .set("cy", cy)
                    .set("r", 2.5);

                active_group = active_group.add(circle);
            }
            document = document.add(active_group);
        }
    } else {
        // No classification - draw all points in red (default)
        let mut sites_group = Group::new()
            .set("fill", "red")
            .set("stroke", "none");

        for center in hole_centers {
            let cx = (center.x - min_x) * scale + margin;
            let cy = (center.y - min_y) * scale + margin;

            let circle = Circle::new()
                .set("cx", cx)
                .set("cy", cy)
                .set("r", 2.0);

            sites_group = sites_group.add(circle);
        }

        document = document.add(sites_group);
    }

    svg::save(output_path, &document)?;

    Ok(())
}