use geo::{Coord, Polygon, Point, Line};
use geo::algorithm::contains::Contains;
use geo::algorithm::euclidean_distance::EuclideanDistance;

/// Generate hexagonal close-packed (HCP) grid points within a polygon
///
/// # Arguments
/// * `polygon` - The boundary polygon to fill with points
/// * `pitch` - The center-to-center spacing between points (in micrometers)
/// * `clearance` - Minimum distance from polygon boundary (in micrometers)
///
/// # Returns
/// Vector of coordinates representing the HCP grid points
pub fn generate_hcp_grid(
    polygon: &Polygon<f64>,
    pitch: f64,
    clearance: f64,
) -> Vec<Coord<f64>> {
    let mut points = Vec::new();

    // Get bounding box
    let coords: Vec<_> = polygon.exterior().coords().collect();
    if coords.is_empty() {
        return points;
    }

    let min_x = coords.iter().map(|c| c.x).fold(f64::INFINITY, f64::min);
    let max_x = coords.iter().map(|c| c.x).fold(f64::NEG_INFINITY, f64::max);
    let min_y = coords.iter().map(|c| c.y).fold(f64::INFINITY, f64::min);
    let max_y = coords.iter().map(|c| c.y).fold(f64::NEG_INFINITY, f64::max);

    println!("\n=== HCP Grid Generation ===");
    println!("Bounding box: ({:.3}, {:.3}) to ({:.3}, {:.3})", min_x, min_y, max_x, max_y);
    println!("Pitch: {:.3} μm", pitch);
    println!("Clearance: {:.3} μm", clearance);

    // HCP grid parameters
    let dx = pitch;
    let dy = pitch * (3.0_f64.sqrt() / 2.0); // Height of equilateral triangle

    // Generate grid
    let mut row = 0;
    let mut y = min_y;

    while y <= max_y {
        let x_offset = if row % 2 == 0 { 0.0 } else { dx / 2.0 };
        let mut x = min_x + x_offset;

        while x <= max_x {
            let coord = Coord { x, y };

            // Check if point is inside polygon
            if polygon.contains(&coord) {
                // Check clearance from boundary if needed
                if clearance > 0.0 {
                    if distance_to_boundary(&coord, polygon) >= clearance {
                        points.push(coord);
                    }
                } else {
                    points.push(coord);
                }
            }

            x += dx;
        }

        y += dy;
        row += 1;
    }

    println!("Generated {} HCP grid points", points.len());

    points
}

/// Compute minimum distance from a point to the polygon boundary
fn distance_to_boundary(point: &Coord<f64>, polygon: &Polygon<f64>) -> f64 {
    let pt = Point::from(*point);
    let exterior = polygon.exterior();

    let mut min_dist = f64::INFINITY;

    // Check distance to each edge of the polygon
    let coords: Vec<_> = exterior.coords().collect();
    for i in 0..coords.len() - 1 {
        let line = Line::new(*coords[i], *coords[i + 1]);
        let dist = pt.euclidean_distance(&line);
        if dist < min_dist {
            min_dist = dist;
        }
    }

    min_dist
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::LineString;

    #[test]
    fn test_hcp_grid_square() {
        // Create a simple square polygon
        let square = Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 10.0),
                (0.0, 10.0),
                (0.0, 0.0),
            ]),
            vec![],
        );

        let points = generate_hcp_grid(&square, 2.0, 0.0);

        // Should have roughly (10/2)^2 * 1.1 points (HCP packs ~10% more than square grid)
        assert!(points.len() > 20);
        assert!(points.len() < 35);

        // All points should be inside the square
        for point in &points {
            assert!(square.contains(point));
        }
    }
}
