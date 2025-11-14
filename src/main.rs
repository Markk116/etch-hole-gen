mod hcp_grid;
mod cvt;
mod dxf_output;
mod svg_output;

use anyhow::{Context, Result};
use dxf::entities::*;
use dxf::Drawing;
use geo::{Coord, LineString, Polygon};
use std::env;

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <input.dxf> [options]", args[0]);
        eprintln!("\nOptions:");
        eprintln!("  --pitch <value>         Hole pitch in μm (default: 2.0)");
        eprintln!("  --diameter <value>      Hole diameter in μm (default: 1.0)");
        eprintln!("  --iterations <value>    Max CVT iterations (default: 50)");
        eprintln!("  --threshold <value>     Convergence threshold (default: 0.001)");
        eprintln!("  --output <path>         Output DXF path (default: output.dxf)");
        eprintln!("  --svg <path>            Optional SVG visualization path");
        std::process::exit(1);
    }

    // Parse command line arguments
    let input_path = &args[1];
    let mut pitch = 2.0;
    let mut hole_diameter = 1.0;
    let mut max_iterations = 50;
    let mut convergence_threshold = 0.001;
    let mut output_path = "output.dxf".to_string();
    let mut svg_path: Option<String> = None;

    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--pitch" => {
                pitch = args[i + 1].parse()?;
                i += 2;
            }
            "--diameter" => {
                hole_diameter = args[i + 1].parse()?;
                i += 2;
            }
            "--iterations" => {
                max_iterations = args[i + 1].parse()?;
                i += 2;
            }
            "--threshold" => {
                convergence_threshold = args[i + 1].parse()?;
                i += 2;
            }
            "--output" => {
                output_path = args[i + 1].clone();
                i += 2;
            }
            "--svg" => {
                svg_path = Some(args[i + 1].clone());
                i += 2;
            }
            _ => {
                eprintln!("Unknown option: {}", args[i]);
                i += 1;
            }
        }
    }

    println!("╔════════════════════════════════════════════════════════╗");
    println!("║   Etch Hole Pattern Generator                          ║");
    println!("║   Centroidal Voronoi Tessellation (CVT)                ║");
    println!("╚════════════════════════════════════════════════════════╝");
    println!("\n📁 Input: {}", input_path);
    println!("📁 Output: {}", output_path);
    if let Some(ref svg) = svg_path {
        println!("📁 SVG: {}", svg);
    }
    println!("\n⚙️  Parameters:");
    println!("   Pitch: {} μm", pitch);
    println!("   Hole diameter: {} μm", hole_diameter);
    println!("   Max iterations: {}", max_iterations);
    println!("   Convergence threshold: {}", convergence_threshold);

    // Step 1: Load DXF and extract polygons
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Step 1: Loading DXF file");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let drawing = Drawing::load_file(input_path)
        .context("Failed to load DXF file")?;

    let polygons = extract_polygons(&drawing)?;

    if polygons.is_empty() {
        anyhow::bail!("No closed polygons found in DXF file!");
    }

    // Use the first polygon as the boundary
    let boundary = &polygons[0];

    println!("✓ Loaded polygon with {} vertices", boundary.exterior().coords().count());

    // Step 2: Generate HCP grid
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Step 2: Generating HCP grid");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let initial_points = hcp_grid::generate_hcp_grid(boundary, pitch, 0.0);

    if initial_points.is_empty() {
        anyhow::bail!("Failed to generate initial HCP grid!");
    }

    println!("✓ Generated {} initial points", initial_points.len());

    // Step 3: Compute CVT
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Step 3: Computing CVT (Lloyd's algorithm)");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let (optimized_points, variance_history) = cvt::compute_cvt(
        initial_points,
        boundary,
        max_iterations,
        convergence_threshold,
    )?;

    println!("✓ CVT optimization complete");
    println!("   Final points: {}", optimized_points.len());
    if !variance_history.is_empty() {
        println!("   Initial variance: {:.6}", variance_history[0]);
        println!("   Final variance: {:.6}", variance_history[variance_history.len() - 1]);
        let improvement = (1.0 - variance_history[variance_history.len() - 1] / variance_history[0]) * 100.0;
        println!("   Improvement: {:.2}%", improvement);
    }

    // Step 4: Write DXF output
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("Step 4: Writing DXF output");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    dxf_output::write_dxf_with_holes(
        &output_path,
        boundary,
        &optimized_points,
        hole_diameter,
    )?;

    println!("✓ DXF file written successfully");

    // Step 5: Optional SVG export
    if let Some(svg_path) = svg_path {
        println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        println!("Step 5: Writing SVG visualization");
        println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

        svg_output::write_svg(
            &svg_path,
            boundary,
            &optimized_points,
            hole_diameter,
            10.0, // 10 pixels per micrometer
        )?;

        println!("✓ SVG file written successfully");
    }

    println!("\n╔════════════════════════════════════════════════════════╗");
    println!("║   ✓ Processing complete!                               ║");
    println!("╚════════════════════════════════════════════════════════╝");

    Ok(())
}

fn extract_polygons(drawing: &Drawing) -> Result<Vec<Polygon<f64>>> {
    let mut polygons = Vec::new();
    let mut open_polylines: Vec<Vec<Coord<f64>>> = Vec::new();

    println!("\n=== Analyzing Entities ===");

    for entity in drawing.entities() {
        match &entity.specific {
            EntityType::LwPolyline(polyline) => {
                let is_closed = (polyline.flags & 1) != 0;
                println!("Found LWPOLYLINE with {} vertices, closed: {}",
                         polyline.vertices.len(), is_closed);

                let coords: Vec<Coord<f64>> = polyline.vertices.iter()
                    .map(|v| Coord { x: v.x, y: v.y })
                    .collect();

                if is_closed && coords.len() >= 3 {
                    let line_string = LineString::new(coords);
                    let polygon = Polygon::new(line_string, vec![]);
                    polygons.push(polygon);
                } else if coords.len() >= 2 {
                    open_polylines.push(coords);
                }
            }
            EntityType::Polyline(polyline) => {
                let is_closed = (polyline.flags & 1) != 0;
                println!("Found POLYLINE, closed: {}", is_closed);
            }
            EntityType::Line(line) => {
                println!("Found LINE from ({:.3}, {:.3}) to ({:.3}, {:.3})",
                         line.p1.x, line.p1.y, line.p2.x, line.p2.y);
            }
            EntityType::Circle(circle) => {
                println!("Found CIRCLE at ({:.3}, {:.3}) with radius {:.3}",
                         circle.center.x, circle.center.y, circle.radius);
            }
            EntityType::Arc(arc) => {
                println!("Found ARC at ({:.3}, {:.3}) with radius {:.3}",
                         arc.center.x, arc.center.y, arc.radius);
            }
            _ => {}
        }
    }

    // Try to assemble open polylines into closed polygons
    if !open_polylines.is_empty() {
        println!("\n=== Open Polyline Segments ===");
        println!("Found {} open polyline segments", open_polylines.len());

        println!("\n=== Attempting to Assemble Open Polylines ===");
        let assembled = assemble_open_polylines(open_polylines)?;
        polygons.extend(assembled);
    }

    if polygons.is_empty() {
        println!("\nWarning: No closed polygons found!");
    }

    Ok(polygons)
}

fn assemble_open_polylines(segments: Vec<Vec<Coord<f64>>>) -> Result<Vec<Polygon<f64>>> {
    if segments.is_empty() {
        return Ok(vec![]);
    }

    let mut all_points: Vec<Coord<f64>> = Vec::new();
    for seg in segments {
        all_points.extend(seg);
    }

    println!("\nCollected {} total points from segments", all_points.len());

    let mut chain = Vec::new();
    let mut used = vec![false; all_points.len()];

    chain.push(all_points[0]);
    used[0] = true;

    while chain.len() < all_points.len() {
        let current = chain.last().unwrap();
        let mut best_idx = None;
        let mut best_dist = f64::INFINITY;

        for (i, point) in all_points.iter().enumerate() {
            if !used[i] {
                let dist = distance(current, point);
                if dist < best_dist {
                    best_dist = dist;
                    best_idx = Some(i);
                }
            }
        }

        if let Some(idx) = best_idx {
            chain.push(all_points[idx]);
            used[idx] = true;
        } else {
            break;
        }
    }

    println!("Built chain with {} points", chain.len());

    let line_string = LineString::new(chain);
    let polygon = Polygon::new(line_string, vec![]);

    Ok(vec![polygon])
}

fn distance(p1: &Coord<f64>, p2: &Coord<f64>) -> f64 {
    ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
}
