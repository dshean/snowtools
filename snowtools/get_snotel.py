#! /usr/bin/env python

"""
Utility to fetch, filter and plot SNOTEL data

Uses ulmo and CUAHSI db
"""

#TODO
#Find date of peak snow depth and SWE
#Add site names to site_list
#Use geom to for spatial filter rather than bbox extent
#Add support for all variable types, generalize handling of scaling and error for each variable (currently hardcoded for SNWD and WTEQ
#Could split this up with multiprocessing rather than running serially - one job per site
#map_plot should update tick labels with mapped coords, not pixel coords
#add subplot below map plot with site elevation vs. x-coord
#Get NED for padded extent
#Shaded relief for map background

import sys
import os
import argparse
from datetime import datetime

import pytz
import numpy as np
import matplotlib.pyplot as plt
import ulmo
from ulmo.util import convert_datetime
from osgeo import gdal

from pygeotools.lib import malib, timelib, geolib, iolib

#URL for query
wsdlurl = "http://worldwater.byu.edu/interactive/snotel/services/index.php/cuahsi_1_1.asmx?WSDL" 

def get_all_lonlat(outdir='.'):
    """
    This will fetch site code, lat, lon for all SNOTEL sites
    Only needs to be run once.
    """
    csv_fn = os.path.join(outdir, 'snotel_lonlat.csv')
    if not os.path.exists(csv_fn):
        print("Generating list of SNOTEL sites and coordinates")
        sites = ulmo.cuahsi.wof.get_sites(wsdlurl)
        lon = []
        lat = []
        code = []
        z = []
        for k,v in sites.iteritems():
            lon.append(float(v['location']['longitude']))
            lat.append(float(v['location']['latitude']))
            code.append(int(v['code']))
            z.append(float(v['elevation_m']))

        out = np.array(zip(code,lon,lat))
        np.savetxt(csv_fn, out, delimiter=',', fmt='%i,%0.5f,%0.5f')
    else:
        out = np.loadtxt(csv_fn, delimiter=',', dtype=None)
    return out

def site_filter_extent(extent, srs=geolib.wgs_srs, pad=None):
    """
    Filter available sites for a given lat/lon extent
    """
    sites = get_all_lonlat()
    sites_srs = geolib.wgs_srs
    if not srs.IsSame(sites_srs):
        print("Converting SNOTEL lat/lon to input ds projection")
        #This returns (x,y,z) coordinate arrays
        sites_proj = np.array(geolib.cT_helper(sites[:,1], sites[:,2], 0, sites_srs, srs)).T
        #Replace the original lon and lat coordinates with projected x and y
        sites[:,1:3] = sites_proj[:,0:2]
    #print(extent)
    #print(sites)
    if pad is not None:
        print("Padding original extent by: %s km" % pad)
        #Convert to meters
        pad *= 1000.
        extent = geolib.pad_extent(extent, width=pad)
        #print(extent)
    valid_idx = ((sites[:,1] > extent[0]) & (sites[:,1] < extent[2]) & (sites[:,2] > extent[1]) & (sites[:,2] < extent[3]))
    valid_sites = sites[valid_idx]
    #Only return site codes, not lat/lon
    #valid_sites = valid_sites[:,0].astype(int)
    if valid_sites.size == 0:
        valid_sites = None
    return valid_sites 

def site_filter_extent_ds(ds, pad=None):
    """
    Filter available sites for a given dataset
    """
    snotel_srs = geolib.wgs_srs
    ds_srs = geolib.get_ds_srs(ds)
    extent = geolib.ds_extent(ds)
    #extent = geolib.ds_extent(ds, snotel_srs)
    #geom = geolib.get_outline(ds)
    return site_filter_extent(extent, ds_srs, pad)
    
def get_series_dt(series, strptime_fmt='%Y-%m-%dT%H:%M:%S'):
    """
    Get datetime series
    """
    #ts = [convert_datetime(vd['date_time_utc']).replace(tzinfo=pytz.utc) for vd in series['values']]
    ts = [datetime.strptime(vd['date_time_utc'], strptime_fmt) for vd in series['values']]
    return np.array(ts, dtype=np.datetime64)

def get_series_val(series):
    """
    Get value series
    """
    # Create a clean timeseries list of (dt,val) tuples
    val = [float(vd['value']) for vd in series['values']]
    val = np.ma.masked_equal(val, -9999)
    val = np.ma.masked_equal(val, 0.0)
    return val

def map_plot(site_list, ds):
    a = iolib.ds_getma(ds)
    clim = malib.calcperc(a, (2,98))
    mX = site_list[:,1]
    mY = site_list[:,2]
    pX, pY = geolib.mapToPixel(mX, mY, ds.GetGeoTransform())
    #f, ax = plt.subplots(1, figsize=(6,6), subplot_kw={'aspect':'equal', 'adjustable':'box-forced'})
    f, ax = plt.subplots(1, figsize=(6,6), subplot_kw={'aspect':'equal'})
    im = ax.imshow(a, vmin=clim[0], vmax=clim[1], cmap='inferno')
    ax.set_facecolor('0.5')
    from imview.lib import pltlib
    pltlib.add_scalebar(ax, geolib.get_res(ds)[0])
    ax.scatter(pX, pY, s=16, facecolors='w', edgecolors='k')
    for i, lbl in enumerate(site_list[:,0]):
        bbox=dict(boxstyle='round,pad=0.1', fc='k', alpha=0.7)
        ax.annotate(str(int(lbl)), xy=(pX[i], pY[i]), xytext=(0, 4), textcoords='offset points', fontsize=8, color='w', bbox=bbox)
    return f

def getparser():
    parser = argparse.ArgumentParser(description="Identify, download, and process SNOTEL records") 
    parser.add_argument('-dt_start', default=None, type=int, help='Start timestamp (format: YYYYMMDD), leave as None for earliest record')
    parser.add_argument('-dt_end', default=None, type=int, help='End timestamp (format: YYYYMMDD), leave as None for latest record')
    parser.add_argument('-outdir', default=os.getcwd(), help='Directory to store intermediate products (default: %(default)s)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-fn', type=str, help='Raster filename to match (e.g., YYYYMMDD_raster.tif)')
    #parser.add_argument('-dt_pad', type=int, default=366, help='Combine data from this many days before and after target date (default: %(default)s)')
    group.add_argument('-extent', default=None, type=float, nargs=4, metavar=('MINLON', 'MINLAT', 'MAXLON', 'MAXLAT'), help='Spatial extent for query')
    group.add_argument('-stack_fn', type=str, help='DEM stack filename to match (e.g., YYYYMMDD_YYYYMMDD_stack_n.tif)')
    parser.add_argument('-extent_pad', type=float, help='Amount to padding for extent, in km')
    #Incremental precip, cumulative precip
    #PRCP (mm), PREC (mm), SNWD (cm), TAVG, TMAX, TMIN, WTEQ (mm)
    vlist_choices = ['PRCP', 'PREC', 'SNWD', 'TAVG', 'TMAX', 'TMIN', 'WTEQ']
    parser.add_argument('-vlist', nargs='+', type=str, choices=vlist_choices, default=['SNWD', 'WTEQ'], help='SNOTEL variables to query')
    return parser

#def main():
parser = getparser()
args = parser.parse_args()

#Set start/end date range for data
#dt_start = datetime(1932,1,1)
dt_start = args.dt_start
#dt_end = datetime.now()
dt_end = args.dt_end
if dt_start is not None:
    dt_start = datetime.strptime(args.dt_start, '%Y%m%d')
if dt_end is not None:
    dt_end = datetime.strptime(args.dt_end, '%Y%m%d')

#Clean this up
if args.fn is not None:
    if os.path.exists(args.fn):
        fn = args.fn
        ds = gdal.Open(fn)
        site_list = site_filter_extent_ds(ds, pad=args.extent_pad) 
elif args.extent is not None:
    site_list = site_filter_extent(extent, pad=args.extent_pad)
elif args.stack_fn is not None:
    #DEM stack, can be used to plot lines on SNOTEL time series
    stack = malib.DEMStack(stack_fn=args.stack_fn)
    dem_dt = stack.date_list
    ds = stack.get_ds()
    site_list = site_filter_extent_ds(ds, pad=args.extent_pad) 
else:
    sys.exit("Must provide valid raster filename or lat/lon extent")

#sitename = 'baker'
#site_list = [999, 909, 1011, 910]
#sitename = 'gm'
#site_list = [622, 682]

if site_list is None:
    sys.exit("No valid sites identified")

vlist = args.vlist

#Accuracy of measurements, in cm
#https://www.wcc.nrcs.usda.gov/snotel/snotel_sensors.html
sigma_factor = 3
snwd_precision = sigma_factor*1.27/100.
wteq_precision = sigma_factor*0.254/100.

print("\nSite codes: %s" % ', '.join(map(str,site_list[:,0].astype(int))))
print("Start date: %s" % dt_start)
print("End date: %s" % dt_end)
print("Variables: %s\n" % ','.join(vlist))

d = {}
for n, site in enumerate(site_list):
    sitecode = int(site[0])
    print('Processing site %i of %i: %i' % ((n+1), len(site_list), sitecode))

    sitekey = 'SNOTEL:%i' % sitecode
    #site = ulmo.cuahsi.wof.get_site_info(wsdlurl, sitekey)

    #Get first variable, use to set dates
    v = vlist[0]
    sitev = 'SNOTEL:%s' % v
    print(sitev)
    series = ulmo.cuahsi.wof.get_values(wsdlurl, sitekey, sitev, start=dt_start, end=dt_end)
    dt = get_series_dt(series)
    d[sitecode] = {'dt':dt}
    val = get_series_val(series)
    d[sitecode][v] = val

    for v in vlist[1:]:
        sitev = 'SNOTEL:%s' % v
        print(sitev)
        series = ulmo.cuahsi.wof.get_values(wsdlurl, sitekey, sitev, start=dt_start, end=dt_end)
        #dt = series['values']['date_time_utc']
        #vals = series['values']['value']
        val = get_series_val(series)
        #Looks like these are not always updated simultaneously, make sure the records are same length
        #Should probably just query both dt and vals simultaneously, rather than assume all variables are same length
        if val.size != dt.size:
            val = val[0:dt.size]
        d[sitecode][v] = val

    #Convert SNWD to m 
    d[sitecode]['SNWD'] /= 100.
    d[sitecode]['WTEQ'] /= 1000.

    #Mask values less than instrument precision
    d[sitecode]['SNWD'] = np.ma.masked_less(d[sitecode]['SNWD'], snwd_precision)
    d[sitecode]['WTEQ'] = np.ma.masked_less(d[sitecode]['WTEQ'], wteq_precision)

    #Calculate density in g/cc
    rho = (d[sitecode]['WTEQ']/d[sitecode]['SNWD'])
    #Mask density values when snow depth is small, helps avoid bogus density values
    depth_thresh = 0.2
    rho[(d[sitecode]['SNWD'] < depth_thresh)] = np.ma.masked
    d[sitecode]['Density'] = rho
    
vlist.append('Density')

print("Plotting")
ts_f, ts_axa = plt.subplots(len(vlist), 1, sharex=True, figsize=(10,7.5))
for sitecode in d.keys():
    #For some reason, can't subtract datetime from np.datetime64
    dt = d[sitecode]['dt'].astype(datetime)
    for n,v in enumerate(vlist):
        vmed = np.ma.median(d[sitecode][v])
        #vmean = np.ma.mean(d[sitecode][[)
        #lbl = '%s: %0.2f' % (sitecode, vmed)
        lbl = str(sitecode)
        p = ts_axa[n].plot(dt, d[sitecode][v], marker='o', ms=1, linestyle='', label=lbl)
        ts_axa[n].set_ylabel(vlist[n])
    ts_axa[n].axhline(vmed, c=p[0].get_color(), linestyle=':', linewidth=0.5)

ts_axa[0].set_ylabel('Snow Depth (m)')
ts_axa[1].set_ylabel('SWE (m w.e.)')
ts_axa[n].set_ylabel('Density (g/cc)')
ts_axa[n].xaxis_date()
ts_axa[n].set_ylim(0,1.0)
ts_axa[n].legend(prop={'size':8})

#Plot lines for DEM timestamps
if args.stack_fn is not None:
    for dt in dem_dt:
        ts_axa[0].axvline(dt, color='k', alpha=0.2)
        ts_axa[2].axvline(dt, color='k', alpha=0.2)
        for sitecode in d.keys():
            #For some reason, can't subtract datetime from np.datetime64
            dt_list = d[sitecode]['dt'].astype(datetime)
            dt_idx = timelib.get_closest_dt_padded_idx(dt.date(), dt_list, pad=3)
            rho_mean = np.mean(d[sitecode]['Density'][dt_idx])
            print(dt, sitecode, rho_mean)

plt.tight_layout()
ts_f.autofmt_xdate()

map_f = map_plot(site_list, ds)

if False:
    fig_fn = '%s_SNOTEL_ts.pdf' % fn 
    plt.savefig(fig_fn, bbox_inches='tight')
    fig_fn = '%s_SNOTEL_ts.png' % fn 
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

    """
    #Limit time series plot to recent years
    ts_axa[n].set_xlim(datetime(2013,8,1), datetime(2016,6,30))
    fig_fn = '%s_SNOTEL_2013-2016.png' % sitename
    f.set_size_inches(4,7.5)
    plt.tight_layout()
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    """

plt.show()
