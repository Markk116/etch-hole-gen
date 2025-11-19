mod hcp_grid;
mod cvt;
mod dxf_output;
mod svg_output;

use anyhow::{Context, Result};
use clap::Parser;
use dxf::entities::*;
use dxf::Drawing;
use geo::{Coord, LineString, Polygon};
use std::time::Instant;

/// Unit conversion constants
const UM_TO_M: f64 = 1e-6;  // micrometers to meters
const MM_TO_M: f64 = 1e-3;  // millimeters to meters

/// Etch hole pattern generator for nanomechanical resonators
///
/// Uses Centroidal Voronoi Tessellation (CVT) to generate optimally
/// distributed hole patterns within arbitrary membrane geometries.
#[derive(Parser, Debug)]
#[command(name = "etch-hole-gen")]
#[command(version, about, long_about = None)]
struct Args {
    /// Input DXF file containing the membrane boundary
    input: String,

    /// Hole pitch in micrometers
    #[arg(short, long, default_value_t = 2.0)]
    pitch: f64,

    /// Hole diameter in micrometers
    #[arg(short, long, default_value_t = 1.0)]
    diameter: f64,

    /// Edge clearance in micrometers (defaults to pitch if not specified)
    #[arg(short, long)]
    clearance: Option<f64>,

    /// Maximum CVT iterations
    #[arg(short, long, default_value_t = 50)]
    iterations: usize,

    /// Convergence threshold in meters
    #[arg(short, long, default_value_t = 1e-9)]
    threshold: f64,

    /// Output DXF file path
    #[arg(short, long, default_value = "output.dxf")]
    output: String,

    /// Optional SVG visualization output path
    #[arg(long)]
    svg: Option<String>,

    /// Optional Voronoi diagram SVG output path
    #[arg(long)]
    voronoi_svg: Option<String>,

    /// Debug mode: output SVG at each iteration with given prefix
    #[arg(long)]
    debug_svg: Option<String>,

    /// Include boundary outline in output DXF
    #[arg(long)]
    include_outline: bool,
}

fn main() -> Result<()> {
    let start_time = Instant::now();
    let args = Args::parse();

    // Convert user input from micrometers to meters (SI base unit)
    let pitch_m = args.pitch * UM_TO_M;
    let diameter_m = args.diameter * UM_TO_M;
    let clearance_m = args.clearance.map(|c| c * UM_TO_M).unwrap_or(pitch_m/2.0);

    // Print header
    println!("Etch Hole Pattern Generator");
    println!("===========================\n");
    println!("Input:  {}", args.input);
    println!("Output: {}", args.output);
    println!("\nParameters:");
    println!("  Pitch:      {:.2} um ({:.3e} m)", args.pitch, pitch_m);
    println!("  Diameter:   {:.2} um ({:.3e} m)", args.diameter, diameter_m);
    println!("  Clearance:  {:.2} um ({:.3e} m)", args.clearance.unwrap_or(args.pitch/2.0), clearance_m);
    println!("  Iterations: {}", args.iterations);
    println!("  Threshold:  {:.3e} m", args.threshold);

    // Step 1: Load DXF and extract boundary
    println!("\n[1/4] Loading DXF file... [{:.2}s]", start_time.elapsed().as_secs_f64());
    let drawing = Drawing::load_file(&args.input)
        .context("Failed to load DXF file")?;

    let boundary = extract_boundary(&drawing)?;
    let vertex_count = boundary.exterior().coords().count();
    println!("      Loaded boundary with {} vertices", vertex_count);

    // Step 2: Generate initial HCP grid
    println!("\n[2/4] Generating HCP grid... [{:.2}s]", start_time.elapsed().as_secs_f64());
    let initial_points = hcp_grid::generate_hcp_grid(&boundary, pitch_m, clearance_m);

    if initial_points.is_empty() {
        anyhow::bail!("No points generated - check pitch and clearance parameters");
    }
    println!("      Generated {} initial points", initial_points.len());

    // Step 3: Optimize with CVT
    println!("\n[3/4] Running CVT optimization... [{:.2}s]", start_time.elapsed().as_secs_f64());
    let (optimized_points, stats) = cvt::compute_cvt(
        initial_points,
        &boundary,
        args.iterations,
        args.threshold,
        args.debug_svg.as_deref(),
        start_time,
    )?;

    println!("      Completed in {} iterations", stats.iterations_run);
    println!("      Final variance: {:.6e}", stats.final_variance);
    println!("      Final elongation: {:.3}", stats.final_elongation);
    if stats.initial_variance > 0.0 {
        let improvement = (1.0 - stats.final_variance / stats.initial_variance) * 100.0;
        println!("      Improvement: {:.1}%", improvement);
    }

    // Step 4: Write outputs
    println!("\n[4/4] Writing output files... [{:.2}s]", start_time.elapsed().as_secs_f64());

    // Write DXF (convert back to mm)
    dxf_output::write_dxf(&args.output, &boundary, &optimized_points, diameter_m, args.include_outline)?;
    println!("      DXF: {}", args.output);

    // Optional SVG outputs
    if let Some(ref svg_path) = args.svg {
        let scale = 1.0 / pitch_m * 10.0; // ~10 pixels per pitch
        svg_output::write_svg(svg_path, &boundary, &optimized_points, diameter_m, scale)?;
        println!("      SVG: {}", svg_path);
    }

    if let Some(ref voronoi_path) = args.voronoi_svg {
        let scale = 1.0 / pitch_m * 10.0;
        svg_output::write_voronoi_svg(voronoi_path, &boundary, &optimized_points, scale)?;
        println!("      Voronoi SVG: {}", voronoi_path);
    }

    println!("\nDone! Generated {} holes. [{:.2}s]", optimized_points.len(), start_time.elapsed().as_secs_f64());

    Ok(())
}

/// Extract the membrane boundary polygon from a DXF drawing
fn extract_boundary(drawing: &Drawing) -> Result<Polygon<f64>> {
    let mut closed_polygons = Vec::new();
    let mut open_segments: Vec<Vec<Coord<f64>>> = Vec::new();

    for entity in drawing.entities() {
        match &entity.specific {
            EntityType::LwPolyline(polyline) => {
                let is_closed = (polyline.flags & 1) != 0;
                let coords: Vec<Coord<f64>> = polyline.vertices.iter()
                    .map(|v| Coord {
                        x: v.x * MM_TO_M,
                        y: v.y * MM_TO_M
                    })
                    .collect();

                if is_closed && coords.len() >= 3 {
                    let polygon = Polygon::new(LineString::new(coords), vec![]);
                    closed_polygons.push(polygon);
                } else if coords.len() >= 2 {
                    open_segments.push(coords);
                }
            }
            EntityType::Spline(spline) => {
                let points = if !spline.fit_points.is_empty() {
                    &spline.fit_points
                } else {
                    &spline.control_points
                };

                if points.len() >= 2 {
                    let coords: Vec<Coord<f64>> = points.iter()
                        .map(|p| Coord {
                            x: p.x * MM_TO_M,
                            y: p.y * MM_TO_M
                        })
                        .collect();
                    open_segments.push(coords);
                }
            }
            _ => {}
        }
    }

    // Return first closed polygon if available
    if let Some(polygon) = closed_polygons.into_iter().next() {
        return Ok(polygon);
    }

    // Otherwise, try to assemble open segments
    if !open_segments.is_empty() {
        let polygon = assemble_segments(open_segments)?;
        return Ok(polygon);
    }

    anyhow::bail!("No valid boundary found in DXF file")
}

/// Assemble open line segments into a closed polygon
fn assemble_segments(mut segments: Vec<Vec<Coord<f64>>>) -> Result<Polygon<f64>> {
    if segments.is_empty() {
        anyhow::bail!("No segments to assemble");
    }

    // Start chain with first segment
    let mut chain = segments.remove(0);

    // Greedily connect nearest segments
    while !segments.is_empty() {
        let chain_start = *chain.first().unwrap();
        let chain_end = *chain.last().unwrap();

        let mut best_idx = None;
        let mut best_dist = f64::INFINITY;
        let mut best_config = 0;

        for (i, seg) in segments.iter().enumerate() {
            let seg_start = *seg.first().unwrap();
            let seg_end = *seg.last().unwrap();

            // Check all 4 possible connections
            let distances = [
                (distance(&chain_end, &seg_start), 1),   // append forward
                (distance(&chain_end, &seg_end), 2),     // append reversed
                (distance(&seg_end, &chain_start), 3),   // prepend forward
                (distance(&seg_start, &chain_start), 4), // prepend reversed
            ];

            for (dist, config) in distances {
                if dist < best_dist {
                    best_dist = dist;
                    best_idx = Some(i);
                    best_config = config;
                }
            }
        }

        if let Some(idx) = best_idx {
            let seg = segments.remove(idx);
            match best_config {
                1 => chain.extend(&seg[1..]),
                2 => {
                    let reversed: Vec<_> = seg.into_iter().rev().collect();
                    chain.extend(&reversed[1..]);
                }
                3 => {
                    let mut new_chain = seg;
                    new_chain.extend(&chain[1..]);
                    chain = new_chain;
                }
                4 => {
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

    Ok(Polygon::new(LineString::new(chain), vec![]))
}

fn distance(p1: &Coord<f64>, p2: &Coord<f64>) -> f64 {
    ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
}
