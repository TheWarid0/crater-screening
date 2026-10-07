#!/usr/bin/env python3
"""
investigate_gawler.py
---------------------
For each cross-matched Gawler candidate, extract BOTH a DEM tile and a magnetic
tile (windowed reads, no full-file load), show them side by side, and save each
as a small GeoTIFF for QGIS. This lets you check whether the topographic
depression and the magnetic anomaly actually coincide at each location -- the
converging-evidence test that separates a real candidate from coincidence.

Usage:
    python investigate_gawler.py /path/to/gawler_dem.vrt /path/to/gawler_TMI.tif

Outputs per candidate:
    <name>_dem.tif / <name>_mag.tif   -> open in QGIS (add hillshade to the DEM)
    <name>_panel.png                  -> quick side-by-side look

Then in QGIS: hillshade the DEM (look for a circular RIM / central peak, not just a
round low), check the magnetic anomaly is ring/circle shaped and centered on the
depression, and overlay the SARIG geological map to rule out granite plutons /
greenstone belts (igneous bodies are circular AND magnetic -- the main false positive).
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from rasterio.windows import from_bounds, Window
from rasterio.warp import transform as warp_transform
from rasterio.crs import CRS
from scipy.ndimage import uniform_filter

WGS84 = CRS.from_epsg(4326)

# The two cross-matched candidates (circular in BOTH DEM and magnetic)
CANDIDATES = [
    ("cand_1", -33.727, 136.681),   # DEM 19107, MAG 79733
    ("cand_2", -32.996, 134.956),   # DEM 13239, MAG 95053
]

HALF_KM = 40           # half-width of the extracted tile
MAX_TILE_PX = 2000     # decimation cap so tiles stay light


def extract(path, lat, lon, half_km=HALF_KM):
    """Windowed read of a (2*half_km) box around (lat,lon), reprojected into the
    file's CRS, auto-decimated. Returns (data, transform, crs, px_m) or None."""
    with rasterio.open(path) as ds:
        xs, ys = warp_transform(WGS84, ds.crs, [lon], [lat])
        cx, cy = xs[0], ys[0]
        half = half_km / 111.0 if ds.crs.is_geographic else half_km * 1000.0
        L, R, B, T = cx - half, cx + half, cy - half, cy + half
        b = ds.bounds
        if R < b.left or L > b.right or T < b.bottom or B > b.top:
            return None
        win = from_bounds(L, B, R, T, ds.transform)
        full = Window(0, 0, ds.width, ds.height)
        win = win.intersection(full)
        w_px, h_px = int(win.width), int(win.height)
        if w_px < 10 or h_px < 10:
            return None
        dec = max(1, int(max(w_px, h_px) / MAX_TILE_PX))
        data = ds.read(1, window=win,
                       out_shape=(max(1, h_px // dec), max(1, w_px // dec))
                       ).astype("float32")
        tr = ds.window_transform(win) * rasterio.Affine.scale(dec, dec)
        px_m = (ds.res[0] if not ds.crs.is_geographic else ds.res[0] * 111000) * dec
        # mask nodata + implausible fills
        nd = ds.nodata
        valid = np.isfinite(data)
        if nd is not None:
            valid &= (data != nd)
        valid &= (data > -1e6) & (data < 1e7)
        data[~valid] = np.nan
        return data, tr, ds.crs, px_m


def residual(field, px_m, regional_km=30):
    """Remove long-wavelength regional field so the local structure stands out."""
    win = max(5, int(regional_km * 1000 / px_m))
    filled = np.nan_to_num(field, nan=float(np.nanmedian(field)))
    return filled - uniform_filter(filled, size=win, mode="reflect")


def save_geotiff(path, data, transform, crs):
    out = np.where(np.isfinite(data), data, -9999).astype("float32")
    with rasterio.open(path, "w", driver="GTiff", height=out.shape[0],
                       width=out.shape[1], count=1, dtype="float32", crs=crs,
                       transform=transform, nodata=-9999, compress="lzw") as dst:
        dst.write(out, 1)


def main(dem_path, mag_path):
    for name, lat, lon in CANDIDATES:
        print(f"\n=== {name}  ({lat}, {lon}) ===")
        dem = extract(dem_path, lat, lon)
        mag = extract(mag_path, lat, lon)
        if dem is None:
            print("  DEM: outside file bounds"); 
        if mag is None:
            print("  MAG: outside file bounds")
        if dem is None or mag is None:
            continue

        dem_data, dem_tr, dem_crs, dem_px = dem
        mag_data, mag_tr, mag_crs, mag_px = mag
        dem_res = residual(dem_data, dem_px)
        mag_res = residual(mag_data, mag_px)
        print(f"  DEM tile {dem_data.shape} @ ~{dem_px:.0f} m/px")
        print(f"  MAG tile {mag_data.shape} @ ~{mag_px:.0f} m/px")

        # save GeoTIFFs for QGIS
        save_geotiff(f"{name}_dem.tif", dem_data, dem_tr, dem_crs)
        save_geotiff(f"{name}_mag.tif", mag_data, mag_tr, mag_crs)

        # quick-look panel: DEM, DEM residual, MAG residual
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
        ax[0].imshow(dem_data, cmap="terrain"); ax[0].set_title(f"{name} DEM (elevation)")
        ax[1].imshow(dem_res, cmap="RdBu_r"); ax[1].set_title("DEM residual (blue=low)")
        ax[2].imshow(mag_res, cmap="RdBu_r"); ax[2].set_title("MAG residual (anomaly)")
        for a in ax:
            a.axis("off")
            # mark tile center (the candidate point)
            a.plot(a.get_xlim()[1] / 2, a.get_ylim()[0] / 2, "k+", ms=12, mew=1.5)
        plt.tight_layout(); plt.savefig(f"{name}_panel.png", dpi=110); plt.close()
        print(f"  wrote {name}_dem.tif, {name}_mag.tif, {name}_panel.png")

    print("\nNext in QGIS: hillshade each DEM tile (look for a circular RIM / central "
          "peak, not just a round low), check the magnetic anomaly is ring-shaped and "
          "centered on the depression, then overlay the SARIG geology to rule out "
          "granite plutons / greenstone belts.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: python investigate_gawler.py <dem.vrt/tif> <magnetic.tif>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
