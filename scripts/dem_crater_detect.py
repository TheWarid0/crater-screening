#!/usr/bin/env python3
"""
dem_crater_detect.py -- find circular topographic depressions (crater candidates)
in a LARGE Geoscience Australia 5m LiDAR DEM GeoTIFF, WITHOUT loading the whole file.

This is for ELEVATION data (bare-earth DEM), not gravity/magnetic. A crater here is a
circular depression, often with a raised rim and/or central peak. We:
  1. read only a window around each candidate (handles 100s of GB files),
  2. auto-decimate so a tile is never too many pixels to handle,
  3. remove regional slope, detect circular depressions via radial symmetry,
  4. export each tile as a small GeoTIFF for QGIS + a PNG quicklook.

NB: the GA 5m DEM is tiled by UTM/MGA zone and only covers surveyed areas
(coastal zone, Murray-Darling basin, population centres). A candidate must fall
inside this file's zone AND inside an actual survey footprint, or the window is empty.
"""
import sys, numpy as np, matplotlib.pyplot as plt
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform as warp_transform
from rasterio.crs import CRS
from scipy.ndimage import uniform_filter, maximum_filter
from scipy.signal import fftconvolve

WGS84 = CRS.from_epsg(4326)
MAX_TILE_PX = 2500   # cap: decimate reads bigger than this per side

CANDIDATES = [
    ("cand_A", -25.5, 139.5), ("cand_B", -25.5, 115.5), ("cand_C", -33.0, 141.5),
    ("cand_D", -23.0, 123.0), ("cand_E", -26.5, 116.5), ("cand_F", -31.0, 135.0),
]

def inspect(path):
    with rasterio.open(path) as ds:
        px_m = ds.res[0] if not ds.crs.is_geographic else ds.res[0]*111000
        print(f"file: {ds.width}x{ds.height}px | CRS {ds.crs} | {ds.res[0]} units/px (~{px_m:.1f} m)")
        print(f"  bounds {ds.bounds}")
        return ds.crs, px_m

def extract(path, lat, lon, half_km=30):
    with rasterio.open(path) as ds:
        xs, ys = warp_transform(WGS84, ds.crs, [lon], [lat])
        cx, cy = xs[0], ys[0]
        half = half_km/111.0 if ds.crs.is_geographic else half_km*1000.0
        L,R,B,T = cx-half, cx+half, cy-half, cy+half
        b = ds.bounds
        if R < b.left or L > b.right or T < b.bottom or B > b.top:
            return None
        win = from_bounds(L,B,R,T, ds.transform)
        # clip the window to the file's actual pixel grid (no boundless padding)
        full = rasterio.windows.Window(0, 0, ds.width, ds.height)
        win = win.intersection(full)
        w_px, h_px = int(win.width), int(win.height)
        if w_px < 10 or h_px < 10:
            return None
        dec = max(1, int(max(w_px,h_px)/MAX_TILE_PX))
        out_w, out_h = max(1,w_px//dec), max(1,h_px//dec)
        data = ds.read(1, window=win, out_shape=(out_h, out_w)).astype("float32")
        # transform for the decimated tile
        tr = ds.window_transform(win)
        tr = tr * tr.scale(dec, dec)
        px_m = (ds.res[0] if not ds.crs.is_geographic else ds.res[0]*111000) * dec
        # detect "no data": GA DEM uses large negative or NaN outside surveys
        valid = np.isfinite(data) & (data > -1000) & (data < 9000)
        if valid.mean() < 0.05:
            return ("EMPTY",)
        data[~valid] = np.nan
        return (data, tr, ds.crs, px_m)

def detect_depressions(dem, px_m, radii_km=(0.5,1,2,4,8), n_top=3):
    """Radial matched filter for circular DEPRESSIONS (central low vs rim)."""
    dem = np.nan_to_num(dem, nan=float(np.nanmedian(dem)))
    # remove regional slope (large-window mean)
    reg = uniform_filter(dem, size=max(5,int(20000/px_m)), mode="reflect")
    res = dem - reg
    best = np.full(res.shape, -np.inf, "float32")
    rm = res - res.mean()
    for rkm in radii_km:
        r_px = rkm*1000/px_m
        k = int(r_px*2.5)|1
        if k<3 or k>min(res.shape): continue
        yy,xx = np.mgrid[0:k,0:k]; c=k//2
        rr = np.sqrt((xx-c)**2+(yy-c)**2)
        templ = -np.exp(-(rr**2)/(2*(r_px*0.6)**2))   # central low
        templ -= templ.mean()
        resp = fftconvolve(rm, templ[::-1,::-1], mode="same")/(np.linalg.norm(templ)+1e-9)
        best = np.maximum(best, resp)
    sep = max(3, int(2000/px_m))
    mx = (best == maximum_filter(best, size=sep))
    ys,xs = np.where(mx & np.isfinite(best))
    vals = best[ys,xs]; order=np.argsort(vals)[::-1][:n_top]
    return res, [(float(vals[i]), int(xs[i]), int(ys[i])) for i in order]

def save_tile_geotiff(path, data, transform, crs):
    data = np.where(np.isfinite(data), data, -9999).astype("float32")
    with rasterio.open(path,"w",driver="GTiff",height=data.shape[0],width=data.shape[1],
                       count=1,dtype="float32",crs=crs,transform=transform,
                       nodata=-9999,compress="lzw") as dst:
        dst.write(data,1)

def main(path, half_km=30):
    crs, px_m = inspect(path)
    print(f"\nextracting {2*half_km} km tiles (auto-decimated to <= {MAX_TILE_PX}px/side)\n")
    for name, lat, lon in CANDIDATES:
        r = extract(path, lat, lon, half_km=half_km)
        if r is None:
            print(f"{name}: outside this file's zone/bounds"); continue
        if isinstance(r[0], str) and r[0] == "EMPTY":
            print(f"{name}: inside bounds but NO LiDAR survey coverage here (empty)"); continue
        dem, tr, tcrs, tpx = r
        res, peaks = detect_depressions(dem, tpx)
        sc = peaks[0][0] if peaks else float("nan")
        print(f"{name}: tile {dem.shape} @ ~{tpx:.0f}m/px | top depression score={sc:.1f}")
        save_tile_geotiff(f"{name}_dem.tif", dem, tr, tcrs)
        fig,ax=plt.subplots(1,2,figsize=(11,4.5))
        ax[0].imshow(dem,cmap="terrain"); ax[0].set_title(f"{name} DEM")
        ax[1].imshow(res,cmap="RdBu_r"); ax[1].set_title("residual (blue=low)")
        for s,px,py in peaks:
            ax[1].plot(px,py,"o",mfc="none",mec="lime",ms=14,mew=2)
        for a in ax: a.axis("off")
        plt.tight_layout(); plt.savefig(f"{name}_dem.png",dpi=110); plt.close()
    print("\nWrote *_dem.tif (open in QGIS) and *_dem.png (quicklook) for covered candidates.")

if __name__ == "__main__":
    if len(sys.argv)<2:
        print("usage: python dem_crater_detect.py /path/to/DEM_z55.tif [half_km]"); sys.exit(1)
    main(sys.argv[1], half_km=int(sys.argv[2]) if len(sys.argv)>2 else 30)
