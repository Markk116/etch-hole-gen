use anyhow::{Context, Result};
use geo::{Coord, LineString, Polygon};
use usvg::{Tree, Node, NodeKind, TreeParsing};
use std::path::Path;

/// Millimeters to meters conversion (SVG typically uses mm or px)
const MM_TO_M: f64 = 1e-3;

/// Extract boundary polygon from an SVG file
///
/// Returns polygon with coordinates in meters.
pub fn extract_boundary_from_svg(path: &Path) -> Result<Polygon<f64>> {
    let svg_data = std::fs::read(path)
        .context("Failed to read SVG file")?;

    let opt = usvg::Options::default();
    let tree = Tree::from_data(&svg_data, &opt)
        .context("Failed to parse SVG")?;

    let mut polygons = Vec::new();

    // Walk the tree to find paths
    extract_paths_recursive(&tree.root, &mut polygons);

    if polygons.is_empty() {
        anyhow::bail!("No paths found in SVG file");
    }

    // Return the largest polygon (by vertex count for now)
    polygons.into_iter()
        .max_by_key(|p| p.exterior().coords().count())
        .ok_or_else(|| anyhow::anyhow!("No valid polygon found"))
}

fn extract_paths_recursive(node: &Node, polygons: &mut Vec<Polygon<f64>>) {
    for child in node.children() {
        match *child.borrow() {
            NodeKind::Path(ref path) => {
                if let Some(polygon) = path_to_polygon(path) {
                    polygons.push(polygon);
                }
            }
            NodeKind::Group(_) => {
                extract_paths_recursive(&child, polygons);
            }
            _ => {}
        }
    }
}

fn path_to_polygon(path: &usvg::Path) -> Option<Polygon<f64>> {
    let mut coords = Vec::new();

    for segment in path.data.segments() {
        match segment {
            usvg::tiny_skia_path::PathSegment::MoveTo(p) => {
                coords.push(Coord {
                    x: p.x as f64 * MM_TO_M,
                    y: p.y as f64 * MM_TO_M,
                });
            }
            usvg::tiny_skia_path::PathSegment::LineTo(p) => {
                coords.push(Coord {
                    x: p.x as f64 * MM_TO_M,
                    y: p.y as f64 * MM_TO_M,
                });
            }
            usvg::tiny_skia_path::PathSegment::QuadTo(_p1, p) => {
                // Approximate quadratic bezier with line to endpoint
                coords.push(Coord {
                    x: p.x as f64 * MM_TO_M,
                    y: p.y as f64 * MM_TO_M,
                });
            }
            usvg::tiny_skia_path::PathSegment::CubicTo(_, _, p) => {
                // Approximate cubic bezier with line to endpoint
                coords.push(Coord {
                    x: p.x as f64 * MM_TO_M,
                    y: p.y as f64 * MM_TO_M,
                });
            }
            usvg::tiny_skia_path::PathSegment::Close => {
                // Close the path
                if let Some(&first) = coords.first() {
                    if coords.last() != Some(&first) {
                        coords.push(first);
                    }
                }
            }
        }
    }

    if coords.len() >= 3 {
        Some(Polygon::new(LineString::new(coords), vec![]))
    } else {
        None
    }
}
