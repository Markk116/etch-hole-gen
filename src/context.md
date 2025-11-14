# Etch Hole Pattern Generator

## Project Goal

Generate optimized hole patterns for nanomechanical resonator fabrication. The holes enable uniform etching of material beneath suspended membrane structures.

## Application Details

**Device Structure:**
- Square membrane: 2mm × 2mm
- 4 tethers: ~50μm wide each
- Tethers connect membrane to frame

**Hole Specifications:**
- Hole diameter: 1μm
- Target spacing: 1μm (2μm pitch center-to-center)
- Total holes per device: ~1 million

**Critical Requirement:**
The membrane-tether interface is especially important due to uneven geometry. Uniform hole distribution is critical for:
- Even undercutting during etching
- Minimizing stress concentrations
- Uniform material removal

## Algorithm Choice: Centroidal Voronoi Tessellation (CVT)

**Goal:** Minimize variance in Voronoi cell sizes across the hole pattern.

**Why CVT:**
- Produces extremely uniform spatial distributions
- Each hole is equidistant from its neighbors (in the Voronoi sense)
- Naturally respects boundary constraints
- Low aspect ratio cells

**Method: Lloyd's Algorithm**
1. Seed with hexagonal close-packed (HCP) grid at ~2μm pitch
2. Iteratively:
   - Compute Voronoi diagram clipped to boundary polygon
   - Calculate centroid of each Voronoi cell
   - Move each point to its cell's centroid
3. Continue until convergence (variance stabilizes)

**Performance Considerations:**
- ~1M points per device
- Expected: 10-100 iterations for convergence
- Each iteration: O(N log N) Voronoi + O(N·M) polygon clipping
- Single-threaded runtime: minutes to hours (acceptable for fab tool)
- Parallelization possible within iterations

## I/O Format

**Input:** DXF file containing closed polygon outline(s) of membrane + tethers

**Output:** DXF file with:
- Original outline
- Hole pattern (circles at optimized positions)
- Ideally using DXF BLOCK/INSERT for instancing (one hole definition, many instances)

## Current Implementation Status

### ✅ Completed
1. **DXF Parser** - Successfully loads DXF files
2. **Polygon Extraction** - Handles broken/segmented paths from Inkscape
3. **Greedy Path Assembly** - Reconstructs closed polygons from disconnected segments using nearest-neighbor chaining

### 🚧 Next Steps
1. **HCP Grid Seeding** - Generate initial hexagonal point distribution inside polygon
2. **CVT Implementation** - Lloyd's algorithm with Voronoi + centroid calculation
3. **Convergence Metrics** - Track variance of cell areas
4. **DXF Output** - Write optimized hole pattern to DXF
5. **Visualization** (optional) - SVG export for quick inspection

## Technical Stack

**Language:** Rust (chosen for performance with large point sets)

**Key Dependencies:**
- `dxf` - DXF file I/O
- `geo` / `geo-types` - Geometric primitives and operations
- Need to add: Voronoi tessellation library (e.g., `voronoi`, `delaunator`, or `spade`)

## Notes

- Boundary clearance TBD (likely want hole centers ≥1μm from edges so hole circumference doesn't intersect boundary)
- May need special handling for narrow tethers
- If CVT is too slow, fallback options: Poisson disk sampling or optimized HCP grid