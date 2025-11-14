use anyhow::{Context, Result};
use dxf::entities::*;
use dxf::Drawing;
use geo::{Coord, LineString, Polygon};
use std::env;

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    
    if args.len() < 2 {
        eprintln!("Usage: {} <path_to_dxf>", args[0]);
        std::process::exit(1);
    }

    let dxf_path = &args[1];
    println!("Loading DXF from: {}", dxf_path);
    
    let drawing = Drawing::load_file(dxf_path)
        .context("Failed to load DXF file")?;
    
    println!("\n=== DXF File Info ===");
    println!("Version: {:?}", drawing.header.version);
    println!("Number of entities: {}", drawing.entities().count());
    
    // Extract and display polygons
    let polygons = extract_polygons(&drawing)?;
    
    println!("\n=== Extracted Polygons ===");
    println!("Found {} closed polygon(s)", polygons.len());
    
    for (i, poly) in polygons.iter().enumerate() {
        println!("\nPolygon {}:", i + 1);
        println!("  Exterior points: {}", poly.exterior().points().count());
        println!("  Interior rings: {}", poly.interiors().len());
        
        // Calculate bounding box
        let coords: Vec<_> = poly.exterior().coords().collect();
        if !coords.is_empty() {
            let min_x = coords.iter().map(|c| c.x).fold(f64::INFINITY, f64::min);
            let max_x = coords.iter().map(|c| c.x).fold(f64::NEG_INFINITY, f64::max);
            let min_y = coords.iter().map(|c| c.y).fold(f64::INFINITY, f64::min);
            let max_y = coords.iter().map(|c| c.y).fold(f64::NEG_INFINITY, f64::max);
            
            println!("  Bounding box: ({:.6}, {:.6}) to ({:.6}, {:.6})", min_x, min_y, max_x, max_y);
            println!("  Width: {:.6}, Height: {:.6}", max_x - min_x, max_y - min_y);
        }
        
        // Display first few points
        println!("  First 5 points:");
        for (j, coord) in poly.exterior().coords().take(5).enumerate() {
            println!("    {}: ({:.6}, {:.6})", j, coord.x, coord.y);
        }
    }
    
    Ok(())
}

fn extract_polygons(drawing: &Drawing) -> Result<Vec<Polygon<f64>>> {
    let mut polygons = Vec::new();
    let mut open_polylines: Vec<Vec<Coord<f64>>> = Vec::new();
    
    println!("\n=== Analyzing Entities ===");
    
    for entity in drawing.entities() {
        match &entity.specific {
            EntityType::LwPolyline(polyline) => {
                let is_closed = (polyline.flags & 1) != 0; // Bit 0 indicates closed
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
                let is_closed = (polyline.flags & 1) != 0; // Bit 0 indicates closed
                println!("Found POLYLINE, closed: {}", is_closed);
                
                // For polylines, we need to look for associated VERTEX entities
                // This is more complex and may require iterating through entities differently
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
            _ => {
                // Uncomment to see all entity types:
                // println!("Found other entity type: {:?}", entity.specific);
            }
        }
    }
    
    // Try to assemble open polylines into closed polygons
    if !open_polylines.is_empty() {
        println!("\n=== Open Polyline Segments ===");
        println!("Found {} open polyline segments", open_polylines.len());
        
        for (i, seg) in open_polylines.iter().enumerate() {
            let start = seg.first().unwrap();
            let end = seg.last().unwrap();
            println!("Segment {}: {} points, Start: ({:.6}, {:.6}), End: ({:.6}, {:.6})",
                     i, seg.len(), start.x, start.y, end.x, end.y);
        }
        
        println!("\n=== Attempting to Assemble Open Polylines ===");
        let assembled = assemble_open_polylines(open_polylines)?;
        polygons.extend(assembled);
    }
    
    if polygons.is_empty() {
        println!("\nWarning: No closed polygons found!");
        println!("The DXF may contain separate LINE segments or open polylines.");
        println!("Consider using a CAD tool to ensure the path is a closed LWPOLYLINE.");
    }
    
    Ok(polygons)
}

fn assemble_open_polylines(segments: Vec<Vec<Coord<f64>>>) -> Result<Vec<Polygon<f64>>> {
    if segments.is_empty() {
        return Ok(vec![]);
    }
    
    // Collect all points from all segments in order
    let mut all_points: Vec<Coord<f64>> = Vec::new();
    for seg in segments {
        all_points.extend(seg);
    }
    
    println!("\nCollected {} total points from segments", all_points.len());
    
    // Build a chain by greedily connecting nearest points
    let mut chain = Vec::new();
    let mut used = vec![false; all_points.len()];
    
    // Start with the first point
    chain.push(all_points[0]);
    used[0] = true;
    
    // Keep adding the nearest unused point to the end of the chain
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
            
            if chain.len() % 5 == 0 {
                println!("  Added {} points, last connection distance: {:.6}", chain.len(), best_dist);
            }
        } else {
            break;
        }
    }
    
    println!("Built chain with {} points", chain.len());
    
    // Check closure
    let start = chain.first().unwrap();
    let end = chain.last().unwrap();
    let closure_distance = distance(start, end);
    
    println!("Chain closure distance: {:.6}", closure_distance);
    println!("  Start: ({:.6}, {:.6})", start.x, start.y);
    println!("  End:   ({:.6}, {:.6})", end.x, end.y);
    
    // Create polygon
    let line_string = LineString::new(chain);
    let polygon = Polygon::new(line_string, vec![]);
    
    Ok(vec![polygon])
}

fn distance(p1: &Coord<f64>, p2: &Coord<f64>) -> f64 {
    ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
}

fn points_match(p1: &Coord<f64>, p2: &Coord<f64>, epsilon: f64) -> bool {
    (p1.x - p2.x).abs() < epsilon && (p1.y - p2.y).abs() < epsilon
}