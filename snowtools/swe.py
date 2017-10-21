#! /usr/bin/env python

"""
Convert snow depth to SWE

Can input two DEMs or a DEM and existing snow depth map
"""

import os
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

from imview.lib import pltlib
from pygeotools.lib import warplib, geolib, iolib, malib, timelib, filtlib

#TODO:
#Allow plotting of snotel sites over SWE map - use get_snotel functions
#Move plots to separate functions
#Better timestamp extraction for dem2 vs. dz filenames
#Implement snow models and better density handling
#Add panel for modscag
#Add panel for TOA reflectance from imagery

def getparser():
    parser = argparse.ArgumentParser(description="Identify, download, and process SNOTEL records") 
    parser.add_argument('-outdir', default=os.getcwd(), help='Output directory (default: %(default)s)')
    parser.add_argument('-dem1_fn', type=str, help='DEM(t1) filename')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-dem2_fn', type=str, default=None, help='DEM(t2) filename')
    group.add_argument('-dz_fn', type=str, default=None, help='Snow depth map filename')
    density_choices = ['sturm', 'snotel', 'constant']
    parser.add_argument('-density', nargs=1, type=float, default=None, help='Specify density (g/cc)')
    parser.add_argument('-filter', action='store_true', help='Filter SWE map to remove blunders, smooth, and fill gaps')
    parser.add_argument('-prism', action='store_true', help='Include PRISM precip in analysis/plots')
    return parser

#def main():
parser = getparser()
args = parser.parse_args()

dem1_fn = args.dem1_fn
dem1_ts = timelib.fn_getdatetime(dem1_fn)
res = 'min'
save = True 
#For testing
#res = 64

if args.dem2_fn is not None:
    dem2_fn = args.dem2_fn
    print("Warping DEMs to same res/extent/proj")
    #This will check input param for validity, could do beforehand
    dem1_ds, dem2_ds = warplib.memwarp_multi_fn([dem1_fn, dem2_fn], extent='intersection', res=res, t_srs='first')
    print("Loading input DEMs into masked arrays")
    dem1 = iolib.ds_getma(dem1_ds, 1)
    dem2 = iolib.ds_getma(dem2_ds, 1)
    dem2_ts = timelib.fn_getdatetime(dem2_fn)
    dz = dem2 - dem1
    outprefix = os.path.splitext(os.path.split(dem1_fn)[1])[0]+'_'+os.path.splitext(os.path.split(dem2_fn)[1])[0]
elif args.dz_fn is not None:
    dz_fn = args.dz_fn
    dem1_ds, dz_ds = warplib.memwarp_multi_fn([dem1_fn, dz_fn], extent='intersection', res=res, t_srs='first')
    print("Loading input DEM and Snow depth into masked arrays")
    dem1 = iolib.ds_getma(dem1_ds, 1)
    dz = iolib.ds_getma(dz_ds, 1)
    #Try to pull out second timestamp from dz_fn
    dem2_ts = timelib.fn_getdatetime_list(dz_fn)[-1]
    outprefix = os.path.splitext(os.path.split(dz_fn)[1])[0]

outprefix = os.path.join(args.outdir, outprefix)

#Calculate water year
wy = dem1_ts.year + 1
if dem1_ts.month >= 10:
    wy = dem1_ts.year

#These need to be updated in geolib to use gdaldem API
hs = geolib.gdaldem_mem_ds(dem1_ds, processing='hillshade', returnma=True)
hs_clim = (1,255)

dem_clim = malib.calcperc(dem1, (1,99))
res = geolib.get_res(dem1_ds)[0]

if args.density is None:
    #Attempt to extract from nearby SNOTEL sites for dem_ts
    #Attempt to use model
    #Last resort, use constant value
    rho_s = 0.5
    #rho_s = 0.4
    #rho_s = 0.36

#Convert snow depth to swe
swe = dz * rho_s

if args.filter:
    print("Filtering SWE map")
    #Median filter to remove artifacts
    swe_f = filtlib.rolling_fltr(swe, size=5)
    #Gaussian filter to smooth over gaps
    swe_f = filtlib.gauss_fltr_astropy(swe, size=9)
    swe = swe_f

swe_clim = list(malib.calcperc(swe, (1,99)))
swe_clim[0] = 0
swe_clim = (0, 8)

prism = None
nax = 2
figsize = (8, 4)
if args.prism:
    #This is PRISM 30-year normal winter PRECIP
    prism_fn = '/Users/dshean/data/PRISM_ppt_30yr_normal_800mM2_10-05_winter_cum.tif'
    if os.path.exists(prism_fn):
        prism_ds = warplib.memwarp_multi_fn([prism_fn,], extent=dem1_ds, res=dem1_ds, t_srs=dem1_ds, r='cubicspline')[0]
        #Values are mm, convert to meters
        prism = iolib.ds_getma(prism_ds)/1000.
        #Apply SWE mask, so we are only considering valid pixels
        prism = np.ma.array(prism, mask=np.ma.getmaskarray(swe))
        nax = 3
        figsize = (12, 4)

if True:
    #Map plots
    f, axa = plt.subplots(1, nax, figsize=figsize, sharex=True, sharey=True, subplot_kw={'aspect':'equal', 'adjustable':'box-forced'})
    hs_im = axa[0].imshow(hs, vmin=hs_clim[0], vmax=hs_clim[1], cmap='gray')
    dem_im = axa[0].imshow(dem1, vmin=dem_clim[0], vmax=dem_clim[1], cmap='cpt_rainbow', alpha=0.5)
    axa[0].set_facecolor('k')
    swe_im = axa[1].imshow(swe, vmin=swe_clim[0], vmax=swe_clim[1], cmap='inferno')
    axa[1].set_facecolor('0.3')
    axa[0].set_title('Late Summer %i' % dem1_ts.year, fontdict={'fontsize':8})
    pltlib.add_cbar(axa[0], dem_im, label='Elevation (m WGS84)')
    pltlib.add_scalebar(axa[0], res)
    axa[1].set_title('WY%i (Summer %i to Spring %i Elev. Diff.)' % (wy, dem1_ts.year, (dem1_ts.year+1)), fontdict={'fontsize':8})
    pltlib.add_cbar(axa[1], swe_im, label=r'SWE Estimate (m w.e., $\rho_s$=0.5)')
    if prism is not None:
        #Should use f.add_subplot() here
        prism_im = axa[2].imshow(prism, vmin=swe_clim[0], vmax=swe_clim[1], cmap='inferno')
        axa[2].set_facecolor('0.3')
        axa[2].set_title('~30-year PRISM Normal: Oct-May Precip', fontdict={'fontsize':8})
        pltlib.add_cbar(axa[2], swe_im, label='Cumulative Precip (m w.e.)')
    for ax in axa:
        pltlib.hide_ticks(ax)
    plt.tight_layout()
    if save:
        fig_fn = '%s_WY%i_SWE_maps.png' % (outprefix, wy)
        if prism is not None:
            fig_fn = os.path.splitext(fig_fn)[0]+'_prism.png'
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

#Regression analysis, 2D Histogram plots
#Needs to be cleaned up and tested
if True:
    #Compare SWE with slope and aspect, each should be optional and handled individually
    slope = geolib.gdaldem_mem_ds(dem1_ds, processing='slope', returnma=True)
    aspect = geolib.gdaldem_mem_ds(dem1_ds, processing='aspect', returnma=True)

    #Mask to preserve one set of valid pixels from all datasets
    mask = malib.common_mask([dem1, swe, slope, aspect])

    dem1 = np.ma.array(dem1, mask=mask).compressed()
    #dem_clim = malib.calcperc(dem1, (5,99.9))
    swe = np.ma.array(swe, mask=mask).compressed()
    slope = np.ma.array(slope, mask=mask).compressed()
    slope_clim = malib.calcperc(slope, (0,99))
    aspect = np.ma.array(aspect, mask=mask).compressed()
    aspect_clim = (0., 360.)
    if prism is not None:
        prism = np.ma.array(prism, mask=mask).compressed()
        #prism_clim = malib.calcperc(prism, (0,99))
        prism_clim = swe_clim

    from imview.lib import pltlib
    f, axa = plt.subplots(1, nax+1, figsize=(10,2.5))
    pltlib.plot_2dhist(axa[0], dem1, swe, dem_clim, swe_clim, maxline=False, trendline=False)
    axa[0].set_ylabel('SWE (m w.e.)')
    axa[0].set_xlabel('Elevation (m WGS84)')
    pltlib.plot_2dhist(axa[1], slope, swe, slope_clim, swe_clim, maxline=False, trendline=False)
    axa[1].set_ylabel('SWE (m w.e.)')
    axa[1].set_xlabel('Slope (deg)')
    pltlib.plot_2dhist(axa[2], aspect, swe, aspect_clim, swe_clim, maxline=False, trendline=False)
    axa[2].set_ylabel('SWE (m w.e.)')
    axa[2].set_xlabel('Aspect (deg)')
    if prism is not None:
        #pltlib.plot_2dhist(axa[3], prism, swe, prism_clim, swe_clim, maxline=False, trendline=False)
        #axa[3].set_ylabel('SWE (m w.e.)')
        #axa[3].set_xlabel('PRISM Precip (m w.e.)')
        pltlib.plot_2dhist(axa[3], dem1, prism, dem_clim, prism_clim, maxline=False, trendline=False)
        axa[3].set_ylabel('PRISM Precip (m w.e.)')
        axa[3].set_xlabel('Elevation (m WGS84)')
    for ax in axa:
        ax.set_facecolor('0.3')
    plt.tight_layout()
    if save:
        fig_fn = '%s_WY%i_SWE_hist.png' % (outprefix, wy)
        if prism is not None:
            fig_fn = os.path.splitext(fig_fn)[0]+'_prism.png'
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

plt.show()
