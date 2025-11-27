# Etch Hole Pattern Generator

Generate optimized hole patterns for nanomechanical resonator fabrication. Takes a DXF or SVG outline and fills it with uniformly distributed etch holes.

## Quick Start

```bash
cargo build --release

# Basic usage (2μm pitch, 1μm diameter holes)
./target/release/etch-hole-gen input.dxf -o output.dxf

# With visualization
./target/release/etch-hole-gen input.dxf -o output.dxf --svg preview.svg

# Custom parameters
./target/release/etch-hole-gen input.dxf \
    --pitch 2.5 --diameter 1.2 --clearance 0.8 \
    --iterations 100 -o output.dxf
```

## Usage

```
etch-hole-gen [OPTIONS] <INPUT>
```

**Input formats:** DXF, SVG (must contain closed polygon defining boundary)

**Units:** Pitch, diameter, and clearance are specified in micrometers. Input file units default to millimeters (use `--input-file-unit` for DXF files in other units).

**Viewing the output**: To view the output, especially in oasis format, [KLayout](https://www.klayout.de/) is highly recommended. 

## About

Developed at TU Delft for nanomechanical resonator fabrication. 

## Options

### Geometry
- `-p, --pitch <PITCH>` - Hole pitch in μm (default: 2)
- `-d, --diameter <DIAMETER>` - Hole diameter in μm (default: 1)
- `-c, --clearance <CLEARANCE>` - Edge clearance in μm (default: pitch/2)

### Algorithm
- `-i, --iterations <ITERATIONS>` - Maximum optimization iterations (default: 50)
- `-t, --threshold <THRESHOLD>` - Convergence threshold in meters, on motion between iterations (default: 1e-9)

### Input/Output
- `-u, --input-file-unit <UNIT>` - Input file unit for DXF: m, mm, um, nm (default: mm)
- `-o, --output <PATH>` - Output DXF file (default: output.dxf)
- `--include-outline` - Include boundary outline in output file

### Optional Outputs
- `--svg <PATH>` - Generate SVG visualization
- `--voronoi-svg <PATH>` - Generate Voronoi diagram visualization
- `--debug-svg <PREFIX>` - Output SVG at each iteration (e.g., `debug_` → `debug_000.svg`, `debug_001.svg`, ...)
- `--gds <PATH>` - Generate GDS file
- `--oasis <PATH>` - Generate OASIS file

### Advanced
- `--include-iso-region <TYPE>` - Keep exact HCP pattern in central region (circle or square)
- `--iso-region-size <SIZE>` - Isolation region size in mm (radius for circle, side length for square)
- `--override-point-limit` - Allow more than 10 million points

The 10 million point limit is to protect against the case of generating billions of points when the wrong input unit is specified. 

At any point hitting `ctrl+c` once let's the current iteration finish and then output whatever results it got. Repeating it aborts the program. 

## Output Formats

**DXF**: Boundary on `BOUNDARY` layer, holes as circles on `HOLES` layer

**GDS/OASIS**: Boundary on layer 1, holes as circles on layer 2

**SVG**: Visual preview (boundary in black, holes in red)

## Examples

Standard 2mm × 2mm membrane with 2μm pitch:
```bash
./etch-hole-gen device.dxf --pitch 2 --diameter 1 -o device_holes.dxf
```

High-density pattern with GDS output:
```bash
./etch-hole-gen device.dxf --pitch 1.5 --diameter 0.8 \
    --clearance 0.5 --iterations 100 \
    -o device_holes.dxf --gds device_holes.gds
```

With exact HCP in center (no optimization artifacts):
```bash
./etch-hole-gen device.dxf --pitch 2 --diameter 1 \
    --include-iso-region circle --iso-region-size 0.5 \
    -o device_holes.dxf
```

Debug optimization process:
```bash
./etch-hole-gen device.dxf --debug-svg iter_ --iterations 50
# Generates iter_000.svg, iter_001.svg, ..., iter_050.svg
```

## Algorithm

Initializes points on a hexagonal close-packed grid, then optimizes their positions using Lloyd's algorithm (iterative Voronoi relaxation). Convergence typically takes 20-50 iterations.

Specifying a iso region can speed up convergence by as much as 2x on the right problem. Can sometimes cause some instability, YMMV.

Written in Rust, therefore :fire: blazingly fast :fire:

## Known Limitations

**Polygons with holes**: The tool can handle outlines with interior holes, but treats them as a single polygon with a connecting edge. This creates an irregularity in the hole pattern along this phantom boundary.

**Large point counts**: For geometries that would generate >10M points, use `--override-point-limit` (expect long runtime).

## Technical Details

**Dependencies**: `dxf`, `geo`, `geo-types`, `voronator`, `svg`

**Performance**: ~120,000 points optimize in ~19 seconds for 11 iterations, with iso region ~16 seconds for 15 iterations. 

**Metrics**: Tracks Voronoi cell area variance (measures uniformity) and cell elongation (measures isotropy). Typical improvement: 40-70% variance reduction.

## Roadmap

Currently v1 is done and usable in production. I'd be open to expanding it with different seed grid types, different cutout shapes than holes, etc. But only if somebody actually needs it. So please open an issue if you do.

## License

MIT

Developed by Mark Kalsbeek at [TU Delft PME](https://www.tudelft.nl/me/over/afdelingen/precision-and-microsystems-engineering-pme).
