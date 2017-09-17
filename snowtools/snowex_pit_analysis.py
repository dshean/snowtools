#! /usr/bin/env

"""
Quick snow depth analysis for SnowEx'17 sites
"""

import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from pygeotools.lib import timelib, malib, iolib, warplib, geolib
from imview.lib import pltlib

#Output from snowex_pit_proc.py
csv_fn = 'snowex_pit_out.csv'

#a = np.genfromtxt(csv_fn, delimiter=',', names=True, dtype=None)
#b = a['datetime', 'x_utm13n', 'y_utm13n', 'depth_m', 'density_kgm3']

a = np.loadtxt(csv_fn, delimiter=',', skiprows=1, dtype=object)
b = a[:,[2,3,4,5,6]]
b[:,0] = timelib.dt2o([datetime.strptime(x, '%Y-%m-%d %H:%M:%S') for x in b[:,0]])
b = np.ma.fix_invalid(b.astype(np.float))

#This should be moved to snowex_pit_proc.py
xlim=(100000, 300000)
b[:,1][(b[:,1] > xlim[1]) | (b[:,1] < xlim[0])] = np.ma.masked
ylim=(4100000, 4500000)
b[:,2][(b[:,2] > ylim[1]) | (b[:,2] < ylim[0])] = np.ma.masked

#Only pass pits with valid x and y coord
b = b[b[:,1:3].count(axis=1) == 2]

#Stereo2SWE preliminary products
dem_fn = '/Users/dshean/Documents/UW/SnowEx/preliminary_mos_20170504/gm_8m-tile-0.tif'
hs_fn = '/Users/dshean/Documents/UW/SnowEx/preliminary_mos_20170504/gm_8m-tile-0_hs_az315.tif'
#snowdepth_fn = '/Users/dshean/Documents/UW/SnowEx/preliminary_snowdepth_20170606/snowdepth_20170201-20170317_mos-tile-0.tif'
snowdepth_fn = '/Users/dshean/Documents/UW/SnowEx/preliminary_snowdepth_20170606/snowdepth_tif/snowdepth_20170201-20170317_mos-tile-0_filt5px.tif'

#Load and clip to common extent
dem_ds, hs_ds, snowdepth_ds = warplib.memwarp_multi_fn([dem_fn, hs_fn, snowdepth_fn], extent='union')
dem = iolib.ds_getma(dem_ds)
hs = iolib.ds_getma(hs_ds)
snowdepth = iolib.ds_getma(snowdepth_ds)

#Pixel coordinates of sample sites
x,y = geolib.mapToPixel(b[:,1], b[:,2], dem_ds.GetGeoTransform())
depth = b[:,3]
rho = b[:,4]

#Sample DEM snow depth
samp = geolib.sample(snowdepth_ds, b[:,1], b[:,2], pad=5)

#Filter to throw out samples with significant roughness over sampled area
samp_perc_thresh = 0.3
samp_mask = (samp[:,1]/samp[:,0]) > samp_perc_thresh
depth_diff = depth - samp[:,0]
depth_diff[samp_mask] = np.ma.masked
idx = np.ma.argsort(np.abs(depth_diff))[::-1]
x = x[idx]
y = y[idx]
depth = depth[idx]
rho = rho[idx]
samp = samp[idx]
depth_diff = depth_diff[idx]

#Scatter point size
s = 16
max_depth = 2.5

#Make plots

#Snow depth scatterplot 
f, ax = plt.subplots()
ax.set_aspect('equal')
clim = malib.calcperc(samp[:,1], (5,95))
#ax.plot(depth, samp[:,0], marker='o', ls='', color='k')
ax.scatter(depth, samp[:,0], c=samp[:,1], cmap='gray', vmin=clim[0], vmax=clim[1], s=36)
ax.set_xlabel('Pit snow depth (m)')
ax.set_ylabel('Stereo2SWE snow depth (m)')
ax.plot([0,3],[0,3],color='r', linewidth=0.5)
ax.set_xlim(0,max_depth)
ax.set_ylim(0,max_depth)
fig_fn = 'snowex_gm_snowdepth_scatter.pdf'
f.savefig(fig_fn, dpi=300, bbox_inches='tight')

#Snow depth histogram of signed diff
f, ax = plt.subplots()
ax.hist(depth_diff, bins=64, range=(-max_depth, max_depth), color='k')
malib.print_stats(depth_diff)
ax.set_xlabel('Snow Depth Diff., Pit - DEM (m)')
ax.set_ylabel('Count')
fig_fn = 'snowex_gm_snowdepth_diff_hist.pdf'
f.savefig(fig_fn, dpi=300, bbox_inches='tight')

#Map plot of snow depth
f, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_facecolor('k')
ax.imshow(hs, cmap='gray')
#clim = malib.calcperc(depth, (2,98))
clim = (0,2)
im = ax.imshow(snowdepth, cmap='inferno', clim=clim, alpha=0.7)
ax.set_ylim(dem.shape[0], 0)
ax.set_xlim(0, dem.shape[1])
pltlib.add_cbar(ax, im, label='Snow Depth (m)')
pltlib.add_scalebar(ax, geolib.get_res(dem_ds)[0])
pltlib.hide_ticks(ax)
fig_fn = 'snowex_gm_snowdepth.png'
f.savefig(fig_fn, dpi=300, bbox_inches='tight')
#Now overlay pits
sc = ax.scatter(x, y, s=s, c=depth, cmap='inferno', vmin=clim[0], vmax=clim[1], edgecolors='k')
fig_fn = 'snowex_gm_snowdepth_pitoverlay.png'
f.savefig(fig_fn, dpi=300, bbox_inches='tight')

#Map plot of elevation
f, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_facecolor('k')
ax.imshow(hs, cmap='gray')
clim = malib.calcperc(dem, (2,98))
im = ax.imshow(dem, cmap='cpt_rainbow', clim=clim, alpha=0.5)
ax.set_ylim(dem.shape[0], 0)
ax.set_xlim(0, dem.shape[1])
pltlib.add_cbar(ax, im, label='Elevation (m)')
pltlib.add_scalebar(ax, geolib.get_res(dem_ds)[0])
pltlib.hide_ticks(ax)
fig_fn = 'snowex_gm_elev.png'
f.savefig(fig_fn, dpi=300, bbox_inches='tight')

#Map plot of density 
f, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_facecolor('k')
ax.imshow(hs, cmap='gray')
clim = malib.calcperc(rho, (2,98))
#clim = (200,400)
sc = ax.scatter(x, y, s=s, c=rho, cmap='inferno', vmin=clim[0], vmax=clim[1], edgecolors='k')
ax.set_ylim(dem.shape[0], 0)
ax.set_xlim(0, dem.shape[1])
pltlib.add_cbar(ax, sc, label='Density (kg/m3)')
pltlib.add_scalebar(ax, geolib.get_res(dem_ds)[0])
pltlib.hide_ticks(ax)
fig_fn = 'snowex_gm_rho.png'
f.savefig(fig_fn, dpi=300, bbox_inches='tight')

#Map plot of snow depth difference 
f, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_facecolor('k')
ax.imshow(hs, cmap='gray')
#clim = malib.calcperc(rho, (2,98))
clim = (-1,1)
sc = ax.scatter(x, y, s=s, c=depth_diff, cmap='RdBu', vmin=clim[0], vmax=clim[1], edgecolors='k')
ax.set_ylim(dem.shape[0], 0)
ax.set_xlim(0, dem.shape[1])
pltlib.add_cbar(ax, sc, label='(Pit - Stereo2SWE) snow depth (m)')
pltlib.add_scalebar(ax, geolib.get_res(dem_ds)[0])
pltlib.hide_ticks(ax)
fig_fn = 'snowex_gm_snowdepth_diff.png'
f.savefig(fig_fn, dpi=300, bbox_inches='tight')

plt.show()
