#!/usr/bin/env python3
"""
circular_scan.py -- scan a LARGE geophysical raster (DEM or magnetic grid) for
circular structures (crater candidates), tile by tile, without loading the whole file.

Two modes:
  --mode dem  : finds circular DEPRESSIONS in topography (central low)
  --mode mag  : finds circular ANOMALIES in magnetic data (low OR ring -- impact
                structures often show a central magnetic low / annular pattern)

Shared engine: windowed read -> regional removal -> multi-scale radial matched
filter -> non-max suppression -> ranked candidate list with lat/lon + score.

Designed for the Gawler Craton (SA) but works on any region/raster.
"""
import sys, argparse, numpy as np
import rasterio
from rasterio.windows import Window
from rasterio.warp import transform as warp_transform
from rasterio.crs import CRS
from scipy.ndimage import uniform_filter, maximum_filter
from scipy.signal import fftconvolve

WGS84 = CRS.from_epsg(4326)

# ---- shared engine ---------------------------------------------------------

def regional_residual(field, px_m, regional_km=30):
    win = max(5, int(regional_km*1000/px_m))
    field = np.nan_to_num(field, nan=float(np.nanmedian(field)))
    return field - uniform_filter(field, size=win, mode="reflect")

def matched_filter(res, px_m, radii_km, polarity="low"):
    """Multi-scale radial matched filter. polarity 'low'=central minimum,
    'ring'=annular (high rim, low center)."""
    best = np.full(res.shape, -np.inf, "float32")
    rm = (res - res.mean()).astype("float32")
    for rkm in radii_km:
        r_px = rkm*1000/px_m
        k = int(r_px*2.5)|1
        if k < 3 or k > min(res.shape): continue
        yy, xx = np.mgrid[0:k,0:k]; c=k//2
        rr = np.sqrt((xx-c)**2+(yy-c)**2)
        if polarity == "ring":
            templ = (np.exp(-((rr-r_px)**2)/(2*(r_px*0.25)**2))    # rim high
                     - np.exp(-(rr**2)/(2*(r_px*0.5)**2)))         # center low
        else:
            templ = -np.exp(-(rr**2)/(2*(r_px*0.6)**2))            # central low
        templ -= templ.mean()
        resp = fftconvolve(rm, templ[::-1,::-1], mode="same")/(np.linalg.norm(templ)+1e-9)
        best = np.maximum(best, resp)
    return best

def rank_peaks(score, px_m, min_sep_km, n_top):
    sep = max(3, int(min_sep_km*1000/px_m))
    mx = (score == maximum_filter(score, size=sep)) & np.isfinite(score)
    ys, xs = np.where(mx)
    if len(ys)==0: return []
    vals = score[ys,xs]; order = np.argsort(vals)[::-1][:n_top]
    return [(float(vals[i]), int(xs[i]), int(ys[i])) for i in order]

# ---- tiled scan over the whole raster -------------------------------------

def scan_raster(path, mode, region=None, tile_px=2000, overlap_px=200,
                radii_km=(2,4,8,15,25), n_per_tile=3, n_final=20):
    polarity = "low" if mode=="dem" else "ring"
    with rasterio.open(path) as ds:
        px_m = ds.res[0] if not ds.crs.is_geographic else ds.res[0]*111000
        # determine pixel bounds of the region of interest
        if region:
            lat0,lat1,lon0,lon1 = region
            xs,ys = warp_transform(WGS84, ds.crs, [lon0,lon1], [lat0,lat1])
            cols = sorted([ds.index(xs[0],ys[0])[1], ds.index(xs[1],ys[1])[1]])
            rows = sorted([ds.index(xs[0],ys[0])[0], ds.index(xs[1],ys[1])[0]])
            c0,c1 = max(0,cols[0]), min(ds.width, cols[1])
            r0,r1 = max(0,rows[0]), min(ds.height, rows[1])
        else:
            c0,c1,r0,r1 = 0,ds.width,0,ds.height

        cands=[]
        step = tile_px - overlap_px
        for ry in range(r0, r1, step):
            for rx in range(c0, c1, step):
                w = min(tile_px, c1-rx); h = min(tile_px, r1-ry)
                if w < 50 or h < 50: continue
                win = Window(rx, ry, w, h)
                data = ds.read(1, window=win).astype("float32")
                valid = np.isfinite(data) & (data > -1e4) & (data < 1e5)
                if valid.mean() < 0.2: continue
                data[~valid] = np.nanmedian(data[valid]) if valid.any() else 0.0
                res = regional_residual(data, px_m)
                score = matched_filter(res, px_m, radii_km, polarity)
                for sc, px, py in rank_peaks(score, px_m, min_sep_km=15, n_top=n_per_tile):
                    X, Y = ds.xy(ry+py, rx+px)
                    lon, lat = warp_transform(ds.crs, WGS84, [X], [Y])
                    cands.append((sc, lat[0], lon[0]))
        cands.sort(reverse=True)
        # de-duplicate nearby peaks from overlapping tiles
        out=[]
        for sc,la,lo in cands:
            if all(np.hypot(la-o[1], lo-o[2])>0.2 for o in out):
                out.append((sc,la,lo))
            if len(out)>=n_final: break
        return out, px_m

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("path"); ap.add_argument("--mode", choices=["dem","mag"], required=True)
    ap.add_argument("--region", nargs=4, type=float, metavar=("LAT0","LAT1","LON0","LON1"))
    args = ap.parse_args()
    region = tuple(args.region) if args.region else None
    cands, px_m = scan_raster(args.path, args.mode, region=region)
    print(f"resolution ~{px_m:.0f} m/px | top circular candidates:")
    print(f"{'rank':4s} {'score':>9s}  lat        lon")
    for i,(sc,la,lo) in enumerate(cands,1):
        print(f"{i:4d} {sc:9.1f}  {la:9.4f}  {lo:9.4f}")


def cross_match(dem_cands, mag_cands, tol_deg=0.15, top=15):
    """Find locations flagged by BOTH the DEM and magnetic scans -- the strongest
    candidates (circular in topography AND in magnetics). Returns combined ranking."""
    matches=[]
    for ds_, dla, dlo in dem_cands:
        for ms_, mla, mlo in mag_cands:
            if np.hypot(dla-mla, dlo-mlo) <= tol_deg:
                # combined score: rank-normalized product favors agreement
                matches.append({
                    "lat": (dla+mla)/2, "lon": (dlo+mlo)/2,
                    "dem_score": ds_, "mag_score": ms_,
                    "combined": ds_*ms_})
                break
    matches.sort(key=lambda m: -m["combined"])
    return matches[:top]

if __name__ == "__main__" and False:
    pass
