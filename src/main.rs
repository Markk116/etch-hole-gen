mod hcp_grid;
mod cvt;
mod dxf_output;
mod svg_output;
mod svg_input;
mod gds_output;
mod oasis_output;

use crate::cvt::IsolationRegion;

use anyhow::{Context, Result};
use clap::Parser;
use dxf::entities::*;
use dxf::Drawing;
use geo::{Coord, LineString, Polygon};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

/// Unit conversion constants
const NM_TO_M: f64 = 1e-9;  // nanometers to meters
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

    /// Edge clearance in micrometers (defaults to pitch/2 if not specified)
    #[arg(short, long)]
    clearance: Option<f64>,

    /// Maximum CVT iterations
    #[arg(short, long, default_value_t = 50)]
    iterations: usize,

    /// Convergence threshold in meters
    #[arg(short, long, default_value_t = 1e-9)]
    threshold: f64,

    /// Input file unit, only important for dxf files
    #[arg(short='u', long, default_value = "mm")]
    input_file_unit: String,

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

    /// Include boundary outline in output
    #[arg(long)]
    include_outline: bool,

    /// Include isolation region to keep HCP exact. Can be circle or square.
    #[arg(long)]
    include_iso_region: Option<String>,

    /// Iso region size in mm. Radius for circle type, side length for square type.
    #[arg(long)]
    iso_region_size: Option<f64>,

    /// Optional GDS output file path
    #[arg(long)]
    gds: Option<String>,

    /// Optional OASIS output file path
    #[arg(long)]
    oasis: Option<String>,

    /// Optional override point-limit break of 10 million points
    #[arg(long)]
    override_point_limit: bool,
}

fn main() -> Result<()> {
    let start_time = Instant::now();
    let args = Args::parse();

    // Set up interrupt handler
    let interrupted = Arc::new(AtomicBool::new(false));
    let interrupted_clone = interrupted.clone();
    ctrlc::set_handler(move || {
        if interrupted_clone.load(Ordering::SeqCst) {
            eprintln!("\nSecond interrupt received, exiting immediately!");
            std::process::exit(130); // 128 + SIGINT(2)
        }
        interrupted_clone.store(true, Ordering::SeqCst);
        eprintln!("\nInterrupt received, finishing current iteration... (press Ctrl-C again to force exit)");
    }).expect("Error setting Ctrl-C handler");

    // Convert user input from micrometers to meters (SI base unit)
    let pitch_m = args.pitch * UM_TO_M;
    let diameter_m = args.diameter * UM_TO_M;
    let clearance_m = args.clearance.map(|c| c * UM_TO_M).unwrap_or(pitch_m/2.0);

    // Scale input 
    let input_unit = args.input_file_unit;
    let input_scale = match input_unit.as_str() {
        "m" => 1.,
        "mm" => MM_TO_M,
        "um" => UM_TO_M,
        "nm" => NM_TO_M,
        unit => panic!("Input unit unknown: {unit:?}"),
    };

    let maybe_iso_region: Option<IsolationRegion> = match args.include_iso_region.as_deref() {
        Some("circle") => {
            if let Some(size) = args.iso_region_size {
                Some(IsolationRegion::Circle {
                    center: Coord { x: 0.0, y: 0.0 },  // Will be auto-detected from boundary centroid
                    radius: size * MM_TO_M,
                })
            } else {
                eprintln!("Warning: circle isolation region requested but no size provided");
                None
            }
        }
        Some("square") => {
            if let Some(size) = args.iso_region_size {
                Some(IsolationRegion::Square {
                    center: Coord { x: 0.0, y: 0.0 },  // Will be auto-detected from boundary centroid
                    side_length: size * MM_TO_M,
                })
            } else {
                eprintln!("Warning: square isolation region requested but no size provided");
                None
            }
        }
        Some(unknown) => {
            eprintln!("Warning: unknown isolation region type '{}', ignoring", unknown);
            None
        }
        None => None,
    };

    // Print header
    println!("Etch Hole Pattern Generator");
    println!("===========================\n");
    println!("Input:      {}", args.input);
    println!("Input Unit: {}", &input_unit);
    println!("Output:     {}", args.output);

    println!("\nParameters:");
    println!("  Pitch:      {:.2} um ({:.3e} m)", args.pitch, pitch_m);
    println!("  Diameter:   {:.2} um ({:.3e} m)", args.diameter, diameter_m);
    println!("  Clearance:  {:.2} um ({:.3e} m)", args.clearance.unwrap_or(args.pitch/2.0), clearance_m);
    println!("  Iterations: {}", args.iterations);
    println!("  Threshold:  {:.3e} m", args.threshold);

    // Step 1: Load input file and extract boundary
    let input_path = Path::new(&args.input);
    let extension = input_path.extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase())
        .unwrap_or_default();

    println!("\n[1/4] Loading input file... [{:.2}s]", start_time.elapsed().as_secs_f64());

    let boundary = match extension.as_str() {
        "svg" => {
            svg_input::extract_boundary_from_svg(input_path)?
        }
        "dxf" | _ => {
            let drawing = Drawing::load_file(&args.input)
                .context("Failed to load DXF file")?;
            extract_boundary_from_dxf(&drawing, &input_scale)?
        }
    };

    let vertex_count = boundary.exterior().coords().count();
    println!("      Loaded boundary with {} vertices", vertex_count);
    let bounding_box = cvt::compute_bounding_box(&boundary);
    let diff_x = bounding_box.1.x - bounding_box.0.x;
    let diff_y = bounding_box.1.y - bounding_box.0.y;
    println!("      Bounding box size: {0:.3}mm by {1:.3}mm", diff_x*1000., diff_y*1000.);

    let max_points = (diff_x * diff_y / (pitch_m * pitch_m)).floor();
    println!("      Initial grid may contain up to {} points.", max_points);
    if max_points > 10e6 && !args.override_point_limit {
        panic!("\nThat's too many points. Add flag '--override-point-limit' to do it anyway.")
    }

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
        clearance_m,
        args.iterations,
        args.threshold,
        args.debug_svg.as_deref(),
        start_time,
        &interrupted,
        maybe_iso_region,
        pitch_m
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

    // Optional GDS output
    if let Some(ref gds_path) = args.gds {
        gds_output::write_gds(
            Path::new(gds_path),
            &boundary,
            &optimized_points,
            diameter_m,
            args.include_outline,
        )?;
        println!("      GDS: {}", gds_path);
    }

    // Optional OASIS output
    if let Some(ref oasis_path) = args.oasis {
        oasis_output::write_oasis(
            Path::new(oasis_path),
            &boundary,
            &optimized_points,
            diameter_m,
            args.include_outline,
        )?;
        println!("      OASIS: {}", oasis_path);
    }

    println!("\nDone! Generated {} holes. [{:.2}s]", optimized_points.len(), start_time.elapsed().as_secs_f64());

    Ok(())
}

/// Extract the membrane boundary polygon from a DXF drawing
fn extract_boundary_from_dxf(drawing: &Drawing, scale_factor: &f64) -> Result<Polygon<f64>> {
    let mut closed_polygons = Vec::new();
    let mut open_segments: Vec<Vec<Coord<f64>>> = Vec::new();

    for entity in drawing.entities() {
        match &entity.specific {
            EntityType::LwPolyline(polyline) => {
                let is_closed = (polyline.flags & 1) != 0;
                let coords: Vec<Coord<f64>> = polyline.vertices.iter()
                    .map(|v| Coord {
                        x: v.x * scale_factor,
                        y: v.y * scale_factor
                    })
                    .collect();

                if is_closed && coords.len() >= 3 {
                    let polygon = Polygon::new(LineString::new(coords), vec![]);
                    closed_polygons.push(polygon);
                } else if coords.len() >= 2 {
                    open_segments.push(coords);
                }
            }
            EntityType::Polyline(polyline) => {
                let is_closed = (polyline.flags & 1) != 0;
                let coords: Vec<Coord<f64>> = polyline.vertices()
                    .map(|v| Coord {
                        x: v.location.x * scale_factor,
                        y: v.location.y * scale_factor
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
                            x: p.x * scale_factor,
                            y: p.y * scale_factor
                        })
                        .collect();
                    open_segments.push(coords);
                }
            }
            anything => {print!("Unknown EntityType found: {:?}\n", anything)}
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
