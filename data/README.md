# Data

Raster data is **not** included because of its size (several GB) and third-party licensing. Put downloaded files in `data/raw/`, which is git-ignored.

| Dataset | Used in | Source | Licence |
|---|---|---|---|
| XGM2019E gravity model | notebooks 03, 04 | Downloaded automatically by `pyshtools.datasets.Earth.XGM2019E` (ICGEM) | see ICGEM |
| NGDC-720 V3 magnetic model | notebook 03 | Downloaded automatically by `pyshtools.datasets.Earth.NGDC_720_V3` | NOAA NCEI |
| SRTM 1 arc-second Global DEM | notebook 05, scripts | [USGS EarthExplorer](https://earthexplorer.usgs.gov/). Mosaic the tiles for lat −34…−30, lon 132…137 into `srtm_gawler_mosaic.tif` | public domain |
| AWAGS MAG 2019 TMI grid, 40 m (`Magmap2019-grid-tmi-Cellsize40m-AWAGS_MAG_2019.ers`) | notebook 05, scripts | Geoscience Australia data portal (GADDS) | CC BY 4.0, © Commonwealth of Australia (Geoscience Australia) |

## `australia_craters.csv`

These are the 30 Australian impact structures used as positives. Coordinates are **approximate**: several values are degree-minute figures written as decimals (for example, Henbury −24.34 means 24°34′S). They are good enough for the ~0.25° gravity grids used here, but should be checked against the [Earth Impact Database](http://www.passc.net/EarthImpactDatabase/) before any finer-scale work.
