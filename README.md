# Etch Hole Pattern Generator

Generate optimized hole patterns for nanomechanical resonator fabrication using **Centroidal Voronoi Tessellation (CVT)**.

## Overview

This tool generates uniform etch hole distributions for suspended membrane structures. The holes enable uniform material removal during the undercutting process, which is critical for minimizing stress concentrations and achieving even etching.

**Algorithm:** Lloyd's algorithm for CVT optimization, starting from a hexagonal close-packed (HCP) seed grid.

## Features

- ✅ DXF file I/O for CAD integration
- ✅ Automatic polygon reconstruction from broken paths (Inkscape-compatible)
- ✅ HCP grid initialization
- ✅ CVT optimization using Lloyd's algorithm
- ✅ Convergence metrics and tracking
- ✅ DXF output with boundary and hole layers
- ✅ SVG visualization export

## Installation

```bash
cargo build --release
```

The binary will be available at `target/release/etch-hole-gen`.

## Usage

```bash
./etch-hole-gen <input.dxf> [options]
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--pitch <value>` | Hole pitch (center-to-center spacing) | 2.0 |
| `--diameter <value>` | Hole diameter | 1.0 |
| `--iterations <value>` | Maximum CVT iterations | 50 |
| `--threshold <value>` | Convergence threshold | 0.001 |
| `--output <path>` | Output DXF file path | output.dxf |
| `--svg <path>` | Optional SVG visualization | none |

**Note:** All dimensional units (pitch, diameter) must match the units used in your input DXF file. For example, if your DXF uses millimeters, specify pitch in millimeters.

### Examples

**Standard nanomechanical resonator (2mm membrane, 2μm pitch):**
```bash
# If DXF is in mm, use mm for pitch/diameter
./etch-hole-gen device.dxf --pitch 0.002 --diameter 0.001 \
    --iterations 50 --output device_with_holes.dxf --svg preview.svg
```

**Test with included test frame:**
```bash
./etch-hole-gen test_data/test_frame.dxf \
    --pitch 0.2 --diameter 0.05 --iterations 30 \
    --output output.dxf --svg output.svg
```

## Input Format

**DXF Requirements:**
- File must contain at least one closed polygon (or path segments that can be assembled)
- The polygon defines the membrane + tether outline
- Units should be consistent (typically millimeters)

**Tip:** If exporting from Inkscape, the tool can automatically reconstruct closed polygons from broken path segments.

## Output Format

**DXF Output:**
- Layer `BOUNDARY`: Original membrane/tether outline (LWPOLYLINE)
- Layer `HOLES`: Optimized hole pattern (CIRCLE entities)

**SVG Output (optional):**
- Visual preview with boundary (black) and holes (red, semi-transparent)
- Useful for quick inspection before fabrication

## Algorithm Details

### 1. HCP Grid Seeding
Generates initial hexagonal close-packed point distribution inside the boundary polygon.

### 2. Lloyd's Algorithm (CVT)
Iteratively:
1. Compute Voronoi diagram for current points
2. Calculate centroid of each Voronoi cell
3. Move each point to its cell's centroid
4. Repeat until convergence

**Convergence:** Stops when average point movement falls below threshold or max iterations reached.

**Metrics:** Tracks variance of Voronoi cell areas (lower = more uniform).

### 3. DXF Generation
Writes boundary and optimized hole positions to DXF file with proper layers.

## Performance

- ~1 million points: Minutes to hours (single-threaded)
- ~4,000 points: ~1-2 seconds for 30 iterations
- ~268 points: Sub-second

**Typical convergence:** 20-50 iterations

## Example Results

Test run with `test_frame.dxf` (4mm × 4mm geometry):
- 268 points at 0.2mm pitch
- 30 iterations in <1 second
- **44% improvement** in cell area variance
- Final variance: 0.004838 (vs initial 0.008646)

## Project Structure

```
src/
├── main.rs          # CLI interface and pipeline
├── hcp_grid.rs      # Hexagonal grid generation
├── cvt.rs           # Lloyd's algorithm implementation
├── dxf_output.rs    # DXF file writing
└── svg_output.rs    # SVG visualization
```

## Technical Stack

- **Language:** Rust 2021
- **Key Dependencies:**
  - `dxf` - DXF file I/O
  - `geo` / `geo-types` - Geometric primitives
  - `voronator` - Voronoi diagram computation
  - `svg` - SVG export

## Limitations & Future Work

- [ ] Proper Voronoi cell clipping to boundary (currently simplified)
- [ ] Boundary clearance enforcement (keep holes away from edges)
- [ ] Parallel processing for large point sets
- [ ] Interactive parameter tuning GUI
- [ ] Support for multiple polygons (e.g., with interior holes)

## References

- **Lloyd's Algorithm:** Lloyd, S.P. (1982). "Least squares quantization in PCM"
- **CVT:** Du et al. (1999). "Centroidal Voronoi Tessellations"

## License

MIT

---

**Application:** Nanomechanical resonators, MEMS devices, suspended membrane structures
**Target:** Uniform material undercutting during fabrication
