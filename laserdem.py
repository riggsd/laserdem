#!/usr/bin/env python3
"""
Convert a Digital Elevation Model (DEM) to layered DXF files for laser cutting.

Each output DXF represents one elevation tier, with two layers:
  - CUT: The shape to cut for this tier
  - ETCH: The outline of the next tier (for alignment), if applicable
"""

import argparse
import logging
from pathlib import Path
from typing import Callable, TypeAlias

import ezdxf
import numpy as np
import rasterio
from rasterio import features
from scipy.ndimage import gaussian_filter, zoom
from shapely import Polygon, MultiPolygon, box, transform as shapely_transform
from shapely.geometry import shape
from shapely.ops import unary_union

PolyGeom: TypeAlias = Polygon | MultiPolygon


def load_dem(path: Path) -> tuple[np.ndarray, rasterio.Affine, float | None, rasterio.crs.CRS]:
    """Load DEM and return array, transform, nodata value, and CRS."""
    with rasterio.open(path) as src:
        dem = src.read(1)
        return dem, src.transform, src.nodata, src.crs


def get_valid_mask(dem: np.ndarray, nodata: float | None) -> np.ndarray:
    """Create boolean mask where True = valid elevation data."""
    if nodata is None:
        return np.ones(dem.shape, dtype=bool)
    return ~np.isclose(dem, nodata)


def smooth_dem(dem: np.ndarray, valid_mask: np.ndarray, sigma: float) -> np.ndarray:
    """
    Apply Gaussian smoothing to DEM without blurring across NODATA boundaries.

    Uses normalized convolution: the Gaussian filter is applied to both the
    elevation values and a validity mask (1=valid, 0=NODATA), then we divide.
    This ensures pixels near the boundary are smoothed using only their valid
    neighbors, properly re-weighted.

    Args:
        dem: Elevation array
        valid_mask: Boolean mask where True = valid data
        sigma: Gaussian kernel standard deviation in pixels
    """
    # Set NODATA pixels to 0 (value doesn't matter due to normalization)
    dem_zeroed = np.where(valid_mask, dem, 0.0)
    weights    = valid_mask.astype(float)

    # Smooth both arrays
    dem_smoothed     = gaussian_filter(dem_zeroed, sigma=sigma)
    weights_smoothed = gaussian_filter(weights, sigma=sigma)

    # Divide to get weighted average using only valid neighbors
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(weights_smoothed > 0, dem_smoothed / weights_smoothed, 0.0)

    # Keep original NODATA pixels unchanged
    return np.where(valid_mask, result, dem)


def smooth_dem_optimized(dem: np.ndarray, valid_mask: np.ndarray, sigma: float) -> np.ndarray:
    """
    Optimized Gaussian smoothing using downsample-smooth-upsample strategy.

    For large smoothing radii, downsampling before smoothing is safe because
    the Gaussian filter already eliminates spatial features smaller than ~2-3 x sigma.
    Downsampling only removes detail that smoothing would eliminate anyway.

    Args:
        dem: Elevation array
        valid_mask: Boolean mask where True = valid data
        sigma: Gaussian kernel standard deviation in pixels
    """
    # Determine safe downsample factor: only downsample if sigma justifies it,
    # cap at 4x, and ensure sigma remains >= ~1.7 pixels after downsampling
    downsample_factor = max(1, min(4, int(sigma // 3)))

    if downsample_factor == 1:
        # No benefit from downsampling, use original method
        return smooth_dem(dem, valid_mask, sigma)

    logging.debug(f"    Downsampling by {downsample_factor}x before smoothing")

    # Downsample using block averaging
    h, w = dem.shape
    new_h = h // downsample_factor
    new_w = w // downsample_factor

    # Trim to multiple of downsample_factor
    dem_trimmed  = dem[:new_h * downsample_factor, :new_w * downsample_factor]
    mask_trimmed = valid_mask[:new_h * downsample_factor, :new_w * downsample_factor]

    # Reshape and average (only where valid)
    dem_blocks  = dem_trimmed.reshape(new_h, downsample_factor, new_w, downsample_factor)
    mask_blocks = mask_trimmed.reshape(new_h, downsample_factor, new_w, downsample_factor)

    # Compute weighted average for each block
    mask_float = mask_blocks.astype(float)
    dem_masked = np.where(mask_blocks, dem_blocks, 0.0)

    dem_down     = dem_masked.sum(axis=(1, 3))
    weights_down = mask_float.sum(axis=(1, 3))

    with np.errstate(divide='ignore', invalid='ignore'):
        dem_down = np.where(weights_down > 0, dem_down / weights_down, 0.0)

    mask_down = weights_down > 0

    # Adjust sigma for downsampled resolution
    sigma_adjusted = sigma / downsample_factor

    # Smooth at downsampled resolution
    dem_smoothed = smooth_dem(dem_down, mask_down, sigma_adjusted)

    # Upsample back to original resolution using bilinear interpolation
    # TODO: consider using downsampled raster for feature extraction rather than scaling to original resolution
    zoom_factor   = downsample_factor
    dem_upsampled = zoom(dem_smoothed, zoom_factor, order=1)

    # Handle size mismatch from trimming
    result = dem.copy()
    result[:dem_upsampled.shape[0], :dem_upsampled.shape[1]] = dem_upsampled

    # Preserve original NODATA pixels
    return np.where(valid_mask, result, dem)


def polygonize_mask(mask: np.ndarray, transform: rasterio.Affine) -> PolyGeom:
    """Convert a binary mask to a polygon/multipolygon geometry."""
    mask_uint8 = mask.astype(np.uint8)
    geoms = []
    for geom_dict, value in features.shapes(mask_uint8, transform=transform):
        if value == 1:
            geoms.append(shape(geom_dict))
    if not geoms:
        raise ValueError("No valid geometry found in mask")
    return unary_union(geoms)


def get_tier_elevations(dem: np.ndarray, mask: np.ndarray, interval: float) -> list[float]:
    """
    Compute tier elevations from lowest to highest.
    
    Starts at the interval below the minimum elevation to ensure the first
    tier encompasses the full DEM perimeter.
    """
    valid_elevations = dem[mask]
    min_elev = float(np.min(valid_elevations))
    max_elev = float(np.max(valid_elevations))
    
    # Start one interval below minimum
    start = np.floor(min_elev / interval) * interval
    # Generate elevations up through max
    elevations = []
    e = start
    while e <= max_elev:
        elevations.append(e)
        e += interval
    return elevations


def create_tier_geometry(
    dem: np.ndarray,
    valid_mask: np.ndarray,
    elevation: float,
    transform: rasterio.Affine
) -> PolyGeom | None:
    """Create polygon for all areas >= elevation within valid data area."""
    tier_mask = (dem >= elevation) & valid_mask
    if not np.any(tier_mask):
        return None
    return polygonize_mask(tier_mask, transform)


def build_coordinate_transformer(
    bounds: tuple[float, float, float, float],
    output_max_mm: float
) -> Callable[[np.ndarray], np.ndarray]:
    """
    Build a function that transforms geo coordinates to mm coordinates.
    
    Output is scaled so the largest dimension equals output_max_mm,
    with origin at lower-left corner.
    """
    minx, miny, maxx, maxy = bounds
    width  = maxx - minx
    height = maxy - miny
    scale  = output_max_mm / max(width, height)
    
    def transform_coords(coords: np.ndarray) -> np.ndarray:
        result = coords.copy()
        result[:, 0] = (coords[:, 0] - minx) * scale
        result[:, 1] = (coords[:, 1] - miny) * scale
        return result
    
    return transform_coords


def transform_geometry(
    geom: PolyGeom,
    transformer: Callable[[np.ndarray], np.ndarray]
) -> PolyGeom:
    """Apply coordinate transformation to a geometry."""
    return shapely_transform(geom, transformer)


def simplify_geometry(
    geom: PolyGeom,
    tolerance_mm: float
) -> PolyGeom:
    """Simplify geometry using Douglas-Peucker algorithm."""
    return geom.simplify(tolerance_mm, preserve_topology=True)


def filter_small_polygons(
    geom: PolyGeom,
    min_size_mm: float
) -> PolyGeom | None:
    """
    Remove polygons smaller than a threshold.
    
    Args:
        geom: Input geometry
        min_size_mm: Minimum size as a length; polygons with area < min_size_mm² are removed
    
    Returns:
        Filtered geometry, or None if all polygons were removed
    """
    min_area = min_size_mm ** 2
    
    if isinstance(geom, Polygon):
        return geom if geom.area >= min_area else None
    
    # MultiPolygon: filter constituents
    kept = [p for p in geom.geoms if p.area >= min_area]
    if not kept:
        return None
    if len(kept) == 1:
        return kept[0]
    return MultiPolygon(kept)


def filter_small_holes(
    geom: PolyGeom,
    min_size_mm: float
) -> PolyGeom:
    """
    Remove holes smaller than a threshold.
    
    Args:
        geom: Input geometry
        min_size_mm: Minimum size as a length; holes with area < min_size_mm² are removed
    
    Returns:
        Geometry with small holes removed
    """
    min_area = min_size_mm ** 2
    
    def filter_polygon_holes(poly: Polygon) -> Polygon:
        kept_holes = [
            hole for hole in poly.interiors
            if Polygon(hole).area >= min_area
        ]
        return Polygon(poly.exterior, kept_holes)
    
    if isinstance(geom, Polygon):
        return filter_polygon_holes(geom)
    
    return MultiPolygon([filter_polygon_holes(p) for p in geom.geoms])


def write_ring_to_layer(msp, ring_coords: list, layer_name: str):
    """Write a single closed ring as an LWPOLYLINE."""
    points = list(ring_coords)
    if points and points[0] != points[-1]:
        points.append(points[0])  # Ensure closed
    msp.add_lwpolyline(points, dxfattribs={"layer": layer_name}, close=True)


def write_geometry_to_layer(msp, geom: PolyGeom, layer_name: str):
    """Write polygon(s) to a DXF layer, handling holes and multipolygons."""
    if geom is None or geom.is_empty:
        return
    
    polygons = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
    
    for poly in polygons:
        # Exterior ring
        write_ring_to_layer(msp, poly.exterior.coords, layer_name)
        # Interior rings (holes)
        for interior in poly.interiors:
            write_ring_to_layer(msp, interior.coords, layer_name)


def write_dxf(
    filepath: Path,
    cut_geom: PolyGeom,
    etch_geom: PolyGeom | None,
    bbox_geom: Polygon,
    tier_bbox_geom: Polygon
):
    """Write a DXF file with CUT layer, optional ETCH layer, BBOX1 (full DEM), and BBOX2 (tier) layers."""
    doc = ezdxf.new("R2010")
    doc.units = ezdxf.units.MM
    msp = doc.modelspace()

    doc.layers.add("CUT", color=1)       # Red
    if etch_geom is not None:
        doc.layers.add("ETCH", color=3)  # Green
    doc.layers.add("BBOX1", color=5)     # Blue
    doc.layers.add("BBOX2", color=2)     # Yellow

    write_geometry_to_layer(msp, cut_geom, "CUT")
    if etch_geom is not None:
        write_geometry_to_layer(msp, etch_geom, "ETCH")
    write_geometry_to_layer(msp, bbox_geom, "BBOX1")
    write_geometry_to_layer(msp, tier_bbox_geom, "BBOX2")

    doc.saveas(filepath)


def process_dem(
    dem_path: Path,
    interval: float,
    output_max_mm: float,
    output_dir: Path,
    simplify_tolerance_mm: float = 0.5,
    smooth_radius: float = 0.0,
    min_island_mm: float = 0.0,
    min_hole_mm: float = 0.0
):
    """Main processing pipeline."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load DEM
    logging.info(f"Loading DEM: {dem_path}")
    dem, transform, nodata, crs = load_dem(dem_path)
    valid_mask = get_valid_mask(dem, nodata)

    # Optional smoothing
    if smooth_radius > 0:
        cell_size = abs(transform.a)
        sigma_pixels = smooth_radius / cell_size
        logging.info(f"Smoothing DEM (radius={smooth_radius} map units, σ={sigma_pixels:.1f} pixels)...")
        dem = smooth_dem_optimized(dem, valid_mask, sigma_pixels)

    # Get tier elevations
    elevations = get_tier_elevations(dem, valid_mask, interval)
    logging.info(f"Elevation tiers ({len(elevations)} total): {[float(e) for e in elevations]}")
    
    # Generate tier geometries (in native CRS)
    logging.info("Generating tier geometries...")
    tier_geoms: dict[float, PolyGeom | None] = {}
    for elev in elevations:
        logging.debug(f"  Generating geometry for tier {elev}")
        tier_geoms[elev] = create_tier_geometry(dem, valid_mask, elev, transform)
    
    # Build coordinate transformer from geo coords to mm
    outline     = polygonize_mask(valid_mask, transform)
    transformer = build_coordinate_transformer(outline.bounds, output_max_mm)
    
    # Create bbox geometry in mm coordinates
    bbox_geo = box(*outline.bounds)
    bbox_mm  = transform_geometry(bbox_geo, transformer)
    
    # Transform and simplify all geometries
    logging.info("Transforming and simplifying geometries...")
    tier_geoms_mm: dict[float, PolyGeom | None] = {}
    for elev, geom in tier_geoms.items():
        if geom is not None:
            logging.debug(f"  Processing tier {elev}")
            geom_mm = transform_geometry(geom, transformer)
            geom_mm = simplify_geometry(geom_mm, simplify_tolerance_mm)
            if min_island_mm > 0:
                geom_mm = filter_small_polygons(geom_mm, min_island_mm)
            if geom_mm is not None and min_hole_mm > 0:
                geom_mm = filter_small_holes(geom_mm, min_hole_mm)
            tier_geoms_mm[elev] = geom_mm
        else:
            tier_geoms_mm[elev] = None
    
    # Write DXF files
    logging.info("Writing DXF files...")
    for i, elev in enumerate(elevations):
        cut_geom = tier_geoms_mm[elev]
        if cut_geom is None:
            continue

        # Next tier for etch layer (if not last tier)
        etch_geom = None
        if i + 1 < len(elevations):
            etch_geom = tier_geoms_mm[elevations[i + 1]]

        # Create bounding box for this tier's cut geometry
        tier_bbox_mm = box(*cut_geom.bounds)

        filename = f"tier_{int(elev):04d}.dxf"
        filepath = output_dir / filename
        write_dxf(filepath, cut_geom, etch_geom, bbox_mm, tier_bbox_mm)
        logging.debug(f"  Wrote {filepath.name}")

    logging.info("Done.")


def main():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s\t%(levelname)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logging.getLogger('rasterio').setLevel(logging.WARNING)
    logging.getLogger('ezdxf').setLevel(logging.WARNING)

    parser = argparse.ArgumentParser(description="Convert DEM to layered DXF files for laser cutting.")
    parser.add_argument("dem",         type=Path,  help="Path to input DEM raster")
    parser.add_argument("interval",    type=float, help="Contour interval (in DEM elevation units)")
    parser.add_argument("output_size", type=float, help="Max output dimension in mm")
    parser.add_argument("output_dir",  type=Path,  help="Output directory for DXF files")
    parser.add_argument("--simplify",   type=float, default=0.5, help="Simplification tolerance in mm (default: 0.5)")
    parser.add_argument("--smooth",     type=float, default=0.0, help="Smoothing radius in map units, e.g. meters for UTM (default: 0, no smoothing)")
    parser.add_argument("--min-island", type=float, default=0.0, help="Remove islands smaller than this size in mm (default: 0, keep all)")
    parser.add_argument("--min-hole",   type=float, default=0.0, help="Remove holes smaller than this size in mm (default: 0, keep all)")
    args = parser.parse_args()

    process_dem(args.dem, args.interval, args.output_size, args.output_dir, args.simplify, args.smooth, args.min_island, args.min_hole)


if __name__ == "__main__":
    main()
