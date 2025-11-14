mod hcp_grid;
mod cvt;
mod particle_sim;
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
        eprintln!("  --clearance <value>     Edge clearance in μm (default: pitch)");
        eprintln!("  --method <cvt|particle> Optimization method (default: cvt)");
        eprintln!("  --iterations <value>    Max iterations (default: 50)");
        eprintln!("  --threshold <value>     Convergence threshold (default: 0.001)");
        eprintln!("  --damping <value>       Particle damping factor (default: 0.5)");
        eprintln!("  --output <path>         Output DXF path (default: output.dxf)");
        eprintln!("  --svg <path>            Optional SVG visualization path");
        eprintln!("  --voronoi-svg <path>    Optional Voronoi diagram SVG path");
        std::process::exit(1);
    }

    // Parse command line arguments
    let input_path = &args[1];
    let mut pitch = 2.0;
    let mut hole_diameter = 1.0;
    let mut clearance: Option<f64> = None;  // None means use pitch as default
    let mut method = "cvt".to_string();
    let mut max_iterations = 50;
    let mut convergence_threshold = 0.001;
    let mut damping = 0.5;
    let mut output_path = "output.dxf".to_string();
    let mut svg_path: Option<String> = None;
    let mut voronoi_svg_path: Option<String> = None;

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
            "--clearance" => {
                clearance = Some(args[i + 1].parse()?);
                i += 2;
            }
            "--method" => {
                method = args[i + 1].clone();
                i += 2;
            }
            "--iterations" => {
                max_iterations = args[i + 1].parse()?;
                i += 2;
            }
            "--damping" => {
                damping = args[i + 1].parse()?;
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
            "--voronoi-svg" => {
                voronoi_svg_path = Some(args[i + 1].clone());
                i += 2;
            }
            _ => {
                eprintln!("Unknown option: {}", args[i]);
                i += 1;
            }
        }
    }

    let method_name = match method.as_str() {
        "particle" => "Particle Physics Simulation",
        _ => "Centroidal Voronoi Tessellation (CVT)",
    };

    println!("╔════════════════════════════════════════════════════════╗");
    println!("║   Etch Hole Pattern Generator                          ║");
    println!("║   {:<53}║", method_name);
    println!("╚════════════════════════════════════════════════════════╝");
    println!("\n📁 Input: {}", input_path);
    println!("📁 Output: {}", output_path);
    if let Some(ref svg) = svg_path {
        println!("📁 SVG: {}", svg);
    }
    // Set clearance to pitch if not specified
    let clearance = clearance.unwrap_or(pitch);

    println!("\n⚙️  Parameters:");
    println!("   Method: {}", method);
    println!("   Pitch: {} μm", pitch);
    println!("   Hole diameter: {} μm", hole_diameter);
    println!("   Edge clearance: {} μm", clearance);
    println!("   Max iterations: {}", max_iterations);
    println!("   Convergence threshold: {}", convergence_threshold);
    if method == "particle" {
        println!("   Damping: {}", damping);
    }

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

    let initial_points = hcp_grid::generate_hcp_grid(boundary, pitch, clearance);

    if initial_points.is_empty() {
        anyhow::bail!("Failed to generate initial HCP grid!");
    }

    println!("✓ Generated {} initial points", initial_points.len());

    // Step 3: Compute optimization (CVT or Particle)
    println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    let (optimized_points, _metric_history) = match method.as_str() {
        "particle" => {
            println!("Step 3: Computing particle physics simulation");
            println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

            let (points, energy_history) = particle_sim::compute_particle_distribution(
                initial_points,
                boundary,
                pitch,
                damping,
                max_iterations,
                convergence_threshold,
            )?;

            println!("✓ Particle simulation complete");
            println!("   Final points: {}", points.len());
            if !energy_history.is_empty() {
                println!("   Initial kinetic energy: {:.6}", energy_history[0]);
                println!("   Final kinetic energy: {:.6}", energy_history[energy_history.len() - 1]);
            }

            (points, energy_history)
        }
        _ => {
            println!("Step 3: Computing CVT (Lloyd's algorithm)");
            println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

            let (points, variance_history) = cvt::compute_cvt(
                initial_points,
                boundary,
                max_iterations,
                convergence_threshold,
            )?;

            println!("✓ CVT optimization complete");
            println!("   Final points: {}", points.len());
            if !variance_history.is_empty() {
                println!("   Initial variance: {:.6}", variance_history[0]);
                println!("   Final variance: {:.6}", variance_history[variance_history.len() - 1]);
                let improvement = (1.0 - variance_history[variance_history.len() - 1] / variance_history[0]) * 100.0;
                println!("   Improvement: {:.2}%", improvement);
            }

            (points, variance_history)
        }
    };

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

    // Step 6: Optional Voronoi diagram SVG export
    if let Some(voronoi_svg_path) = voronoi_svg_path {
        println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        println!("Step 6: Writing Voronoi diagram SVG");
        println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

        svg_output::write_voronoi_svg(
            &voronoi_svg_path,
            boundary,
            &optimized_points,
            10.0, // 10 pixels per micrometer
        )?;

        println!("✓ Voronoi SVG file written successfully");
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
            EntityType::Spline(spline) => {
                println!("Found SPLINE with {} fit points, {} control points",
                         spline.fit_points.len(), spline.control_points.len());

                // Use fit points if available, otherwise use control points
                let points = if !spline.fit_points.is_empty() {
                    &spline.fit_points
                } else {
                    &spline.control_points
                };

                if points.len() >= 2 {
                    let coords: Vec<Coord<f64>> = points.iter()
                        .map(|p| Coord { x: p.x, y: p.y })
                        .collect();
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
            other => {
                // Log any other entity types we encounter
                println!("Found other entity type: {:?}", std::mem::discriminant(other));
            }
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

    println!("\n=== Nearest-Neighbor Chain Assembly ===");

    // Track which segments we've used
    let mut remaining_segments: Vec<Vec<Coord<f64>>> = segments;

    let mut polygons = Vec::new();

    while !remaining_segments.is_empty() {
        // Start a new chain with the first remaining segment
        let mut chain = remaining_segments.remove(0);

        println!("Starting new chain with {} points", chain.len());

        // Keep connecting nearest segments until none are left
        while !remaining_segments.is_empty() {
            let chain_start = *chain.first().unwrap();
            let chain_end = *chain.last().unwrap();

            // Find the nearest segment endpoint to either end of our chain
            let mut best_idx = None;
            let mut best_dist = f64::INFINITY;
            let mut best_config = 0; // 1=append_fwd, 2=append_rev, 3=prepend_fwd, 4=prepend_rev

            for (i, seg) in remaining_segments.iter().enumerate() {
                let seg_start = *seg.first().unwrap();
                let seg_end = *seg.last().unwrap();

                // Check all 4 possible connections and pick the closest
                let d1 = distance(&chain_end, &seg_start); // append forward
                let d2 = distance(&chain_end, &seg_end);   // append reversed
                let d3 = distance(&seg_end, &chain_start); // prepend forward
                let d4 = distance(&seg_start, &chain_start); // prepend reversed

                if d1 < best_dist {
                    best_dist = d1;
                    best_idx = Some(i);
                    best_config = 1;
                }
                if d2 < best_dist {
                    best_dist = d2;
                    best_idx = Some(i);
                    best_config = 2;
                }
                if d3 < best_dist {
                    best_dist = d3;
                    best_idx = Some(i);
                    best_config = 3;
                }
                if d4 < best_dist {
                    best_dist = d4;
                    best_idx = Some(i);
                    best_config = 4;
                }
            }

            // Connect the best match
            if let Some(idx) = best_idx {
                let seg = remaining_segments.remove(idx);

                match best_config {
                    1 => {
                        // Append forward
                        println!("  Connecting segment {} to end (forward), distance: {:.6}", idx, best_dist);
                        chain.extend(&seg[1..]);
                    }
                    2 => {
                        // Append reversed
                        println!("  Connecting segment {} to end (reversed), distance: {:.6}", idx, best_dist);
                        let reversed: Vec<_> = seg.into_iter().rev().collect();
                        chain.extend(&reversed[1..]);
                    }
                    3 => {
                        // Prepend forward
                        println!("  Connecting segment {} to start (forward), distance: {:.6}", idx, best_dist);
                        let mut new_chain = seg;
                        new_chain.extend(&chain[1..]);
                        chain = new_chain;
                    }
                    4 => {
                        // Prepend reversed
                        println!("  Connecting segment {} to start (reversed), distance: {:.6}", idx, best_dist);
                        let mut reversed: Vec<_> = seg.into_iter().rev().collect();
                        reversed.extend(&chain[1..]);
                        chain = reversed;
                    }
                    _ => unreachable!(),
                }
            } else {
                break;
            }
        }

        println!("Chain complete with {} points", chain.len());

        // Check if chain is closed
        let chain_start = *chain.first().unwrap();
        let chain_end = *chain.last().unwrap();
        let closure_dist = distance(&chain_start, &chain_end);

        println!("  Start: ({:.6}, {:.6})", chain_start.x, chain_start.y);
        println!("  End:   ({:.6}, {:.6})", chain_end.x, chain_end.y);
        println!("  Closure distance: {:.6}", closure_dist);

        if closure_dist < 1e-3 {
            println!("  ✓ Chain is closed");
        } else {
            println!("  ⚠ Chain has gap (this may cause issues)");
        }

        let line_string = LineString::new(chain);
        let polygon = Polygon::new(line_string, vec![]);
        polygons.push(polygon);
    }

    println!("Assembled {} polygon(s)", polygons.len());

    Ok(polygons)
}

fn distance(p1: &Coord<f64>, p2: &Coord<f64>) -> f64 {
    ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
}
