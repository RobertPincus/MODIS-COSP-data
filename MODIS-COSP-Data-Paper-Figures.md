---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.0
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Dependencies

```{code-cell} ipython3
import datetime
import re 
import copy

import xarray as xr
import numpy as np
import pandas as pd

import seaborn as sns
import colorcet as cc
import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
```

# Graphical definitions

```{code-cell} ipython3
#
# Graphical choices
#
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams["legend.frameon"] = False
mpl.rcParams['axes.labelsize']= 14

saveFigs = True

map_projection = ccrs.Mollweide
```

```{code-cell} ipython3
def global_mean(ds):
    ''' Area-weighted global-mean value of a single field '''
    weights = np.cos(np.abs(ds.latitude * np.arccos(-1)/180.))
    weights = weights/weights.mean()
    return((ds * weights).mean(dim=["latitude", "longitude"]))
```

```{code-cell} ipython3
figDir  = pathlib.Path("../")
dataDir = pathlib.Path("/Users/robert/Codes/MODIS-COSP-Data")
sample_month = "2021-07-01"
```

# Pixel counts - one day, one month

These will require access to/download two example files.

+++

# Cloud mask and retrieval fractions

```{code-cell} ipython3
month = xr.open_dataset(dataDir.joinpath("modis-cosp-scalars.nc"), engine="netcdf4").sel(date=sample_month)
```

```{code-cell} ipython3
sns.set_context("paper")
fig = plt.figure(figsize = (12, 9))
axes = fig.subplots(nrows=2, subplot_kw={'projection': map_projection()})

# Cloud mask fraction  
pl = month.Cloud_Mask_Fraction.plot(ax = axes[0], 
                                    transform=ccrs.PlateCarree(),
                                    cmap = cc.m_gray,
                                    add_colorbar=False)
axes[0].set_title("Cloud mask fraction, July 2021")

pl = month.Cloud_Retrieval_Fraction_Total.plot(ax = axes[1], 
                                    transform=ccrs.PlateCarree(),
                                    cmap = cc.m_gray,
                                    add_colorbar=False)
axes[1].set_title("Cloud retrieval fraction, July 2021")
fig.colorbar(pl, ax=axes.ravel().tolist(), shrink = 0.75, aspect = 15, label="Cloud fraction")

if saveFigs: 
    # fig.savefig(figDir.joinpath("Cloud-fractions.pdf"), dpi=600, transparent=True, bbox_inches = "tight")
    fig.savefig(figDir.joinpath("Cloud-fractions.png"), dpi=150, transparent=True, bbox_inches = "tight")
```

# Differences between mask and retrieval fractions 

Still need to print out mean value of difference

```{code-cell} ipython3
sns.set_context("paper")
fig = plt.figure(figsize = (12, 9))
axes = fig.subplots(nrows=2, subplot_kw={'projection': map_projection()})
cm = cc.m_CET_L19.copy()
cm.set_bad("0.85")

cf_diff = (month.Cloud_Mask_Fraction - month.Cloud_Retrieval_Fraction_Total)
pl = cf_diff.plot(ax = axes[0], transform=ccrs.PlateCarree(),
          vmin = 0, 
          vmax = .4, 
          cmap = cm, 
          add_colorbar=False)
axes[0].set_title("Cloud mask fraction - cloud retrieval fraction, July 2021")
print("Global mean mask fraction - retrieval fraction is {:4.1f}%".format(
      100 * global_mean(cf_diff)))



cf_diff = (month.Cloud_Mask_Fraction - month.Cloud_Retrieval_Fraction_Total - month.Cloud_Retrieval_Fraction_PCL_Total)
pl = cf_diff.plot(ax = axes[1], transform=ccrs.PlateCarree(),
          vmin = 0, 
          vmax = .4, 
          cmap = cm, 
          add_colorbar=False)
axes[1].set_title("Cloud mask fraction - (cloud retrieval + PCL) fraction, July 2021")
print("Global mean mask fraction - (retrieval + PCL) fraction is {:4.1f}%".format(
      100 * global_mean(cf_diff)))

fig.colorbar(pl, ax=axes.ravel().tolist(), shrink = 0.75, aspect = 15, label="Cloud fraction difference")


if saveFigs: 
    # fig.savefig(figDir.joinpath("Cloud-fraction-diff-map.pdf"), dpi=600, transparent=True, bbox_inches = "tight")
    fig.savefig(figDir.joinpath("Cloud-fraction-diff-maps.png"), dpi=150, transparent=True, bbox_inches = "tight")
```

# Hi/mid/low clouds

```{code-cell} ipython3
sns.set_context("paper")
fig = plt.figure(figsize = (12, 12))
axes = fig.subplots(nrows=3, subplot_kw={'projection': map_projection()})

cm = copy.copy(cc.m_fire)
cm.set_bad("0.85")
# Cloud mask fraction
rat = (month.Cloud_Mask_Fraction_High/month.Cloud_Mask_Fraction)
pl = rat.plot(ax = axes[0], 
          transform=ccrs.PlateCarree(),
          vmin = 0, vmax =  1, 
          cmap = cm, 
          add_colorbar=False)
axes[0].set_title("Fraction of clouds that are high, July 2021")


rat = (month.Cloud_Mask_Fraction_Mid/month.Cloud_Mask_Fraction)
pl = rat.plot(ax = axes[1], 
          transform=ccrs.PlateCarree(),
          vmin = 0, vmax =  1, 
          cmap = cm, 
          add_colorbar=False)
axes[1].set_title("Fraction of clouds that are mid, July 2021")

rat = (month.Cloud_Mask_Fraction_Low/month.Cloud_Mask_Fraction)
pl = rat.plot(ax = axes[2], 
          transform=ccrs.PlateCarree(),
          vmin = 0, vmax =  1, 
          cmap = cm, 
          add_colorbar=False)
axes[2].set_title("Fraction of clouds that are low, July 2021")




fig.colorbar(pl, ax=axes.ravel().tolist(), shrink = 0.75, aspect = 15, label="Cloud fraction")

if saveFigs: 
    # fig.savefig(figDir.joinpath("Cloud-fractions.pdf"), dpi=600, transparent=True, bbox_inches = "tight")
    fig.savefig(figDir.joinpath("Cloud-height-partitioning.png"), dpi=150, transparent=True, bbox_inches = "tight")
```

# Liquid/ice clouds

```{code-cell} ipython3
sns.set_context("paper")
fig = plt.figure(figsize = (12, 9))
axes = fig.subplots(nrows=2, subplot_kw={'projection': map_projection()})

cm = copy.copy(cc.m_CET_CBL1)
cm.set_bad("0.85")
# Cloud mask fraction
rat = (month.Cloud_Retrieval_Fraction_Liquid/month.Cloud_Retrieval_Fraction_Total)
pl = rat.plot(ax = axes[0], 
          transform=ccrs.PlateCarree(),
          vmin = 0, vmax =  1, 
          cmap = cm, 
          add_colorbar=False)
axes[0].set_title("Fraction of clouds that are liquid, July 2021")

rat = (month.Cloud_Retrieval_Fraction_Ice/month.Cloud_Retrieval_Fraction_Total)
pl = rat.plot(ax = axes[1], 
          transform=ccrs.PlateCarree(),
          vmin = 0, vmax =  1, 
          cmap = cm, 
          add_colorbar=False)
axes[1].set_title("Fraction of clouds that are ice, July 2021")

fig.colorbar(pl, ax=axes.ravel().tolist(), shrink = 0.75, aspect = 15, label="Cloud fraction")

if saveFigs: 
    # fig.savefig(figDir.joinpath("Cloud-fractions.pdf"), dpi=600, transparent=True, bbox_inches = "tight")
    fig.savefig(figDir.joinpath("Cloud-phase-partitioning.png"), dpi=150, transparent=True, bbox_inches = "tight")
```

# Optical thickness - geometric and arithmetic means

```{code-cell} ipython3
sns.set_context("paper")
fig = plt.figure(figsize = (12, 9))
axes = fig.subplots(nrows=2, subplot_kw={'projection': map_projection()})

cm = copy.copy(cc.m_CET_L18_r)
cm.set_bad("0.85")
#cm = cc.m_CET_L20

# Cloud optical thickness
data = month["Cloud_Optical_Thickness_Total"]
pl = data.plot(ax = axes[0], 
          transform=ccrs.PlateCarree(),
          vmin=0, vmax=40,
          cmap = cm, 
          add_colorbar=False)
axes[0].set_title("Arithemtic mean optical thickness, July 2021")

data = np.power(10, month["Cloud_Optical_Thickness_Log10_Total"])
pl = data.plot(ax = axes[1], 
          transform=ccrs.PlateCarree(),
          vmin=0, vmax=40,
          cmap = cm, 
          add_colorbar=False)
axes[1].set_title("Geometric mean optical thickness, July 2021")

fig.colorbar(pl, ax=axes.ravel().tolist(), shrink = 0.75, aspect = 15, label="Optical thickness")


if saveFigs: 
    # fig.savefig(figDir.joinpath("Cloud-fractions.pdf"), dpi=600, transparent=True, bbox_inches = "tight")
    fig.savefig(figDir.joinpath("Cloud-optical-thickness.png"), dpi=150, transparent=True, bbox_inches = "tight")
```

# How to compute a global-mean joint histogram

+++

# Then make three joint histogram figures

```{code-cell} ipython3
s = "" 
v = "Total"
hname = f"Cloud_Optical_Thickness{s}_{v}_vs_Cloud_Top_Pressure"
fname = dataDir.joinpath(f"modis-cosp-{hname}.nc")
hname = f"Cloud_Optical_Thickness{s}_{v}_vs_Cloud_Top_Pressure"
fname = dataDir.joinpath(f"modis-cosp-{hname}.nc")
gmh = global_mean(xr.open_dataset(fname).sel(date=sample_month))
# Normalize the color bar to the figure with the biggest cloud fractions 
if v == "Total" and s == "": vmax = np.max(gmh[hname])
    
vmax
```

```{code-cell} ipython3
sns.set_context("paper")
fig = plt.figure(figsize = (7.75, 8.9))

cmap = cc.m_blues
axes = fig.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
#
# Joint histograms vs cloud top pressure: (cloudy, partly-cloudy) x (total, liquid, ice)
# 
for i, s in enumerate(["", "_PCL"]): # Top row/subset is cloudy, bottom row/subset is partly-cloudy (PCL)
    for r, v in enumerate(["Total", "Liquid", "Ice"]):
        hname = f"Cloud_Optical_Thickness{s}_{v}_vs_Cloud_Top_Pressure"
        fname = dataDir.joinpath(f"modis-cosp-{hname}.nc")
        gmh = global_mean(xr.open_dataset(fname).sel(date=sample_month))
        # Normalize the color bar to the figure with the biggest cloud fractions 
        if v == "Total" and s == "": vmax = np.max(gmh[hname])
        #    
        # imshow() seems to show the first array dimension from top to bottom and the second dimension 
        # from left to right, so we transpose the DataArray before plotting 
        #
        pl = axes[r, i].imshow(gmh[hname].T, cmap = cmap, vmin=0, vmax=vmax) 
        #
        # Labeling the axes
        #
        if r==2: 
            # Simplified axis label
            axes[r, i].set_xlabel("cloud optical thickness")
            tau_var  = [           k for k in gmh.coords if "optical_thickness"     in k and "bnds" in k][0]
            tau_bnds = gmh[tau_var] 
            axes[r, i].set_xticks(np.arange(len(tau_bnds))-.5)
            axes[r, i].set_xticklabels(tau_bnds.values)
        #
        # y axis label, tick values on left-most panel only
        #
        if i==0: 
            jnt_var  = [           k for k in gmh.coords if "optical_thickness" not in k and "bnds" in k][0]
            axes[r, i].set_ylabel(jnt_var.replace("bnds", "").replace("_", ""))
            jnt_bnds = gmh[jnt_var] 
            axes[r, i].set_yticks(np.arange(len(jnt_bnds))-.5)
            axes[r, i].set_yticklabels(jnt_bnds.values)



axes[0,0].annotate("Total",  (0,0))
axes[1,0].annotate("Liquid", (0,0))
axes[2,0].annotate("Ice",    (0,0))
axes[0,0].annotate("Fully-cloudy",  (4,0))
axes[0,1].annotate("Partly-cloudy", (4,0))

fig.tight_layout()

fig.colorbar(pl, ax=axes.ravel().tolist(), shrink = 0.5, aspect = 15, label="Cloud fraction")

fig.savefig(figDir.joinpath("Tau-pc-histograms.pdf"), dpi=600, transparent=True, bbox_inches = "tight")
```

```{code-cell} ipython3
sns.set_context("paper")
fig = plt.figure(figsize = (10.5, 6.75))
cmap = cc.m_CET_L17 # New colormap because different axes
axes = fig.subplots(ncols=2, nrows=2, sharex=True)
#
# Joint histograms of optical thickness vs particle size 
# 
for i, s in enumerate(["", "_PCL"]): # Top row/subset is cloudy, bottom row/subset is partly-cloudy (PCL)
    for r, v in enumerate(["Liquid", "Ice"]): 
        hname = f"Cloud_Optical_Thickness_vs_Cloud_Particle_Size{s}_{v}"
        fname = dataDir.joinpath(f"modis-cosp-{hname}.nc")
        gmh = global_mean(xr.open_dataset(fname).sel(date=sample_month))
        if v == "Liquid" and s == "": vmax = np.max(gmh[hname])
        #    
        # imshow() seems to show the first array dimension from top to bottom and the second dimension 
        # from left to right, so we transpose the DataArray before plotting 
        #
        pl = axes[r, i].imshow(gmh[hname].T, cmap, vmin=0, vmax=vmax) # Normalization from the figure above
        #
        # Labeling the axes
        #
        # x-axis labels only in bottom row
        if r==1: 
            # Simplified axis label
            axes[r, i].set_xlabel("cloud optical thickness")
            tau_var  = [           k for k in gmh.coords if "optical_thickness"     in k and "bnds" in k][0]
            tau_bnds = gmh[tau_var] 
            axes[r, i].set_xticks(np.arange(len(tau_bnds))-.5)
            axes[r, i].set_xticklabels(tau_bnds.values)

        #
        # y axis labels - need separately for each column 
        #
        jnt_var  = [           k for k in gmh.coords if "optical_thickness" not in k  and "bnds" in k][0]
        axes[r, i].set_ylabel(jnt_var.replace("bin_boundaries", "").replace("_", " ").replace(" pcl", ""))
        jnt_bnds = gmh[jnt_var] 
        axes[r, i].set_yticks(np.arange(len(jnt_bnds))-.5)
        axes[r, i].set_yticklabels(jnt_bnds.values)
        axes[r, i].invert_yaxis()
axes[0,0].annotate("Liquid",  (0,5))
axes[1,0].annotate("Ice",     (0,5))
axes[0,0].annotate("Fully-cloudy",  (4,5))
axes[0,1].annotate("Partly-cloudy", (4,5))


fig.tight_layout()
fig.colorbar(pl, ax=axes.ravel().tolist(), shrink = 0.75, aspect = 15, label="Cloud fraction")
fig.savefig(figDir.joinpath("Tau-re-histograms.pdf"), dpi=600, transparent=True, bbox_inches = "tight")
```

```{code-cell} ipython3
sns.set_context("paper")
fig = plt.figure(figsize = (10.5, 6.75))
cmap = cc.m_CET_CBTL2_r # New colormap because different axes
axes = fig.subplots(ncols=2, nrows=2)
#
# Joint histograms vs cloud top pressure: (cloudy, partly-cloudy) x (total, liquid, ice)
# 
for i, s in enumerate(["", "_PCL"]): # Top row/subset is cloudy, bottom row/subset is partly-cloudy (PCL)
    for r, v in enumerate(["Liquid", "Ice"]): 
        hname = f"Cloud_Water_Path_vs_Cloud_Particle_Size{s}_{v}"
        fname = dataDir.joinpath(f"modis-cosp-{hname}.nc")
        gmh = global_mean(xr.open_dataset(fname).sel(date=sample_month))
        if v == "Liquid" and s == "": vmax = np.max(gmh[hname])
        #    
        # imshow() seems to show the first array dimension from top to bottom and the second dimension 
        # from left to right, so we transpose the DataArray before plotting 
        #
        pl = axes[r, i].imshow(gmh[hname].T, cmap, vmin=0, vmax=vmax) # Normalization from the figure above
        #
        # Labeling the axes
        #
        tau_var  = [           k for k in gmh.coords if "cloud_particle" not in k and "bnds" in k][0]
        #
        # Simplified axis label
        #
        axes[r, i].set_xlabel(tau_var.replace("bin_boundaries", "").replace("_", " ").replace("pcl ", ""))
        tau_bnds = gmh[tau_var] 
        axes[r, i].set_xticks(np.arange(len(tau_bnds))-.5)
        axes[r, i].set_xticklabels(tau_bnds.values)
        #
        # y axis labels - need separately for each column 
        #
        jnt_var  = [           k for k in gmh.coords if "cloud_particle"     in k  and "bnds" in k][0]
        axes[r, i].set_ylabel(jnt_var.replace("bin_boundaries", "").replace("_", " ").replace(" pcl", ""))
        jnt_bnds = gmh[jnt_var] 
        axes[r, i].set_yticks(np.arange(len(jnt_bnds))-.5)
        axes[r, i].set_yticklabels(jnt_bnds.values)
        axes[r, i].invert_yaxis()

        axes[0,0].annotate("Liquid",  (0,5))
axes[1,0].annotate("Ice",     (0,5))
axes[0,0].annotate("Fully-cloudy",  (4,5))
axes[0,1].annotate("Partly-cloudy", (4,5))


fig.tight_layout()
fig.colorbar(pl, ax=axes.ravel().tolist(), shrink = 0.75, aspect = 15, label="Cloud fraction")
fig.savefig(figDir.joinpath("LWP-re-histograms.pdf"), dpi=600, transparent=True, bbox_inches = "tight")
```

```{code-cell} ipython3

```
