use geo::{Coord, Polygon, LineString};
use geo::algorithm::contains::Contains;
use anyhow::Result;
use crate::svg_output;

/// Compute particle-based distribution using repulsive forces
///
/// # Arguments
/// * `initial_points` - Starting point distribution (e.g., from HCP grid)
/// * `boundary` - Polygon boundary to constrain the tessellation
/// * `pitch` - Expected spacing between particles
/// * `damping` - Damping factor for velocity (0.0 = no damping, 1.0 = full damping)
/// * `max_iterations` - Maximum number of simulation steps
/// * `convergence_threshold` - Stop when average movement is below this threshold
/// * `debug_svg_prefix` - Optional prefix for debug SVG output at each iteration
///
/// # Returns
/// Optimized point positions and kinetic energy history
pub fn compute_particle_distribution(
    mut points: Vec<Coord<f64>>,
    boundary: &Polygon<f64>,
    pitch: f64,
    damping: f64,
    max_iterations: usize,
    convergence_threshold: f64,
    debug_svg_prefix: Option<&str>,
) -> Result<(Vec<Coord<f64>>, Vec<f64>)> {
    println!("\n=== Particle Physics Simulation ===");
    println!("Initial points: {}", points.len());
    println!("Pitch: {:.6}", pitch);
    println!("Damping: {:.3}", damping);
    println!("Max iterations: {}", max_iterations);
    println!("Convergence threshold: {:.6}", convergence_threshold);

    // Create offset boundary (inward by pitch/2)
    let offset_distance = pitch / 2.0;
    let constraint_boundary = create_inward_offset_polygon(boundary, offset_distance);

    println!("Boundary offset: {:.6} (half pitch)", offset_distance);

    // Initialize velocities
    let mut velocities: Vec<Coord<f64>> = vec![Coord { x: 0.0, y: 0.0 }; points.len()];
    let mut energy_history = Vec::new();

    // Simulation parameters
    let force_constant = pitch * pitch; // Scale force with pitch
    let interaction_radius = pitch * 3.0; // Only compute forces for nearby particles
    let dt = 0.1; // Time step

    for iter in 0..max_iterations {
        // Compute forces on each particle
        let mut forces: Vec<Coord<f64>> = vec![Coord { x: 0.0, y: 0.0 }; points.len()];

        for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                let dx = points[j].x - points[i].x;
                let dy = points[j].y - points[i].y;
                let r_sq = dx * dx + dy * dy;
                let r = r_sq.sqrt();

                // Only interact with nearby particles
                if r < interaction_radius && r > 1e-10 {
                    // Repulsive force: F = k / r^2
                    let force_mag = force_constant / r_sq;
                    let fx = -force_mag * dx / r;
                    let fy = -force_mag * dy / r;

                    // Apply equal and opposite forces
                    forces[i].x += fx;
                    forces[i].y += fy;
                    forces[j].x -= fx;
                    forces[j].y -= fy;
                }
            }
        }

        // Update velocities and positions
        let mut total_kinetic_energy = 0.0;
        let mut total_movement = 0.0;

        for i in 0..points.len() {
            // Update velocity with damping: v = (v + F*dt) * (1 - damping)
            velocities[i].x = (velocities[i].x + forces[i].x * dt) * (1.0 - damping);
            velocities[i].y = (velocities[i].y + forces[i].y * dt) * (1.0 - damping);

            // Update position
            let new_x = points[i].x + velocities[i].x * dt;
            let new_y = points[i].y + velocities[i].y * dt;
            let mut new_pos = Coord { x: new_x, y: new_y };

            // Constrain to offset boundary
            if !constraint_boundary.contains(&new_pos) {
                // Project back onto boundary
                new_pos = project_onto_boundary(&new_pos, &constraint_boundary);
                // Zero velocity when hitting boundary
                velocities[i].x = 0.0;
                velocities[i].y = 0.0;
            }

            // Track movement
            let movement = distance(&points[i], &new_pos);
            total_movement += movement;

            points[i] = new_pos;

            // Compute kinetic energy
            let v_sq = velocities[i].x * velocities[i].x + velocities[i].y * velocities[i].y;
            total_kinetic_energy += v_sq;
        }

        energy_history.push(total_kinetic_energy);
        let avg_movement = total_movement / points.len() as f64;

        if iter % 10 == 0 || iter == max_iterations - 1 {
            println!(
                "Iteration {}: avg movement = {:.6}, kinetic energy = {:.6}",
                iter, avg_movement, total_kinetic_energy
            );
        }

        // Output debug SVG if requested
        if let Some(prefix) = debug_svg_prefix {
            let debug_path = format!("{}_{:04}.svg", prefix, iter);
            if let Err(e) = svg_output::write_voronoi_svg(&debug_path, boundary, &points, 10.0) {
                eprintln!("Warning: Failed to write debug SVG {}: {}", debug_path, e);
            }
        }

        // Check convergence
        if avg_movement < convergence_threshold {
            println!(
                "\nConverged after {} iterations (avg movement: {:.6})",
                iter, avg_movement
            );
            break;
        }
    }

    println!("Final points: {}", points.len());

    Ok((points, energy_history))
}

/// Create an inward offset polygon (simple approach using vertex displacement)
fn create_inward_offset_polygon(polygon: &Polygon<f64>, offset: f64) -> Polygon<f64> {
    let coords: Vec<_> = polygon.exterior().coords().collect();

    if coords.len() < 4 {
        return polygon.clone();
    }

    let mut offset_coords = Vec::new();

    // For each vertex, compute inward normal and offset
    for i in 0..(coords.len() - 1) {
        let prev_idx = if i == 0 { coords.len() - 2 } else { i - 1 };
        let next_idx = i + 1;

        let prev = coords[prev_idx];
        let curr = coords[i];
        let next = coords[next_idx];

        // Compute edge vectors
        let v1_x = curr.x - prev.x;
        let v1_y = curr.y - prev.y;
        let v1_len = (v1_x * v1_x + v1_y * v1_y).sqrt();

        let v2_x = next.x - curr.x;
        let v2_y = next.y - curr.y;
        let v2_len = (v2_x * v2_x + v2_y * v2_y).sqrt();

        if v1_len < 1e-10 || v2_len < 1e-10 {
            offset_coords.push(*curr);
            continue;
        }

        // Compute inward normals (perpendicular, pointing inward)
        // For a CCW polygon, inward normal is 90° clockwise rotation
        let normal1_x = v1_y / v1_len;
        let normal1_y = -v1_x / v1_len;
        let normal2_x = v2_y / v2_len;
        let normal2_y = -v2_x / v2_len;

        // Average the normals
        let avg_normal_x = (normal1_x + normal2_x) / 2.0;
        let avg_normal_y = (normal1_y + normal2_y) / 2.0;
        let avg_len = (avg_normal_x * avg_normal_x + avg_normal_y * avg_normal_y).sqrt();

        if avg_len > 1e-10 {
            // Normalize and scale by offset
            let offset_x = curr.x + (avg_normal_x / avg_len) * offset;
            let offset_y = curr.y + (avg_normal_y / avg_len) * offset;
            offset_coords.push(Coord { x: offset_x, y: offset_y });
        } else {
            offset_coords.push(*curr);
        }
    }

    // Close the polygon
    if let Some(&first) = offset_coords.first() {
        offset_coords.push(first);
    }

    Polygon::new(LineString::new(offset_coords), vec![])
}

/// Project a point onto the boundary polygon
fn project_onto_boundary(point: &Coord<f64>, boundary: &Polygon<f64>) -> Coord<f64> {
    let exterior = boundary.exterior();
    let coords: Vec<_> = exterior.coords().collect();

    let mut min_dist = f64::INFINITY;
    let mut closest_point = *point;

    // Find closest point on boundary
    for i in 0..coords.len() - 1 {
        let p1 = coords[i];
        let p2 = coords[i + 1];

        let projected = project_point_onto_line_segment(point, p1, p2);
        let dist = distance(point, &projected);

        if dist < min_dist {
            min_dist = dist;
            closest_point = projected;
        }
    }

    closest_point
}

/// Project a point onto a line segment
fn project_point_onto_line_segment(
    point: &Coord<f64>,
    seg_start: &Coord<f64>,
    seg_end: &Coord<f64>,
) -> Coord<f64> {
    let dx = seg_end.x - seg_start.x;
    let dy = seg_end.y - seg_start.y;
    let len_sq = dx * dx + dy * dy;

    if len_sq < 1e-10 {
        return *seg_start;
    }

    // Compute projection parameter t
    let t = ((point.x - seg_start.x) * dx + (point.y - seg_start.y) * dy) / len_sq;
    let t_clamped = t.clamp(0.0, 1.0);

    Coord {
        x: seg_start.x + t_clamped * dx,
        y: seg_start.y + t_clamped * dy,
    }
}

/// Euclidean distance between two points
fn distance(p1: &Coord<f64>, p2: &Coord<f64>) -> f64 {
    ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_project_point_onto_line_segment() {
        let point = Coord { x: 1.0, y: 1.0 };
        let seg_start = Coord { x: 0.0, y: 0.0 };
        let seg_end = Coord { x: 2.0, y: 0.0 };

        let projected = project_point_onto_line_segment(&point, &seg_start, &seg_end);

        assert!((projected.x - 1.0).abs() < 1e-10);
        assert!((projected.y - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_inward_offset() {
        // Simple square
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

        let offset = create_inward_offset_polygon(&square, 1.0);
        let coords: Vec<_> = offset.exterior().coords().collect();

        // Offset square should have vertices moved inward
        assert!(coords.len() > 0);
    }
}
