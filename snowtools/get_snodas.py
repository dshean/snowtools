#! /usr/bin/env python
"""
Utility to identify, download and process SNODAS model snow depth to match an input raster

"""

import sys
import os
import argparse

from osgeo import gdal
import numpy as np

from datetime import datetime, timedelta

from pygeotools.lib import iolib
from pygeotools.lib import warplib
from pygeotools.lib import timelib

def get_snodas(dt, outdir=None, code=1036):
    """Function to fetch and process SNODAS snow depth products for input datetime

    http://nsidc.org/data/docs/noaa/g02158_snodas_snow_cover_model/index.html

    Product codes:
    1036 is snow depth
    1034 is SWE

    filename format: us_ssmv11036tS__T0001TTNATS2015042205HP001.Hdr

    """
    import tarfile
    import gzip
    snodas_ds = None
    snodas_url_str = None
    #Note: unmasked products (beyond CONUS) are only available from 2010-present
    if dt >= datetime(2003,9,30) and dt < datetime(2010,1,1):
        snodas_url_str = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/%Y/%m_%b/SNODAS_%Y%m%d.tar'
        tar_subfn_str_fmt = 'us_ssmv1%itS__T0001TTNATS%%Y%%m%%d05HP001.%s.gz'
    elif dt >= datetime(2010,1,1):
        snodas_url_str = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/unmasked/%Y/%m_%b/SNODAS_unmasked_%Y%m%d.tar'
        tar_subfn_str_fmt = './zz_ssmv1%itS__T0001TTNATS%%Y%%m%%d05HP001.%s.gz'
    else:
        print("No SNODAS data available for input date")

    if snodas_url_str is not None:
        snodas_url = dt.strftime(snodas_url_str)
        snodas_tar_fn = iolib.getfile(snodas_url, outdir=outdir)
        print("Unpacking")
        tar = tarfile.open(snodas_tar_fn)
        #gunzip to extract both dat and Hdr files, tar.gz
        for ext in ('dat', 'Hdr'):
            tar_subfn_str = tar_subfn_str_fmt % (code, ext)
            tar_subfn_gz = dt.strftime(tar_subfn_str)
            tar_subfn = os.path.splitext(tar_subfn_gz)[0]
            print(tar_subfn)
            if outdir is not None:
                tar_subfn = os.path.join(outdir, tar_subfn)
            if not os.path.exists(tar_subfn):
                #Should be able to do this without writing intermediate gz to disk
                tar.extract(tar_subfn_gz)
                with gzip.open(tar_subfn_gz, 'rb') as f:
                    outf = open(tar_subfn, 'wb')
                    outf.write(f.read())
                    outf.close()
                os.remove(tar_subfn_gz)

        #Need to delete 'Created by module comment' line from Hdr, can contain too many characters
        bad_str = 'Created by module comment'
        snodas_fn = tar_subfn
        f = open(snodas_fn)
        output = []
        for line in f:
            if not bad_str in line:
                output.append(line)
        f.close()
        f = open(snodas_fn, 'w')
        f.writelines(output)
        f.close()

        #Return GDAL dataset for extracted product
        #Note: this is entire CONUS
        snodas_ds = gdal.Open(snodas_fn)
    return snodas_ds

def getparser():
    parser = argparse.ArgumentParser(description="Identify, download, and process SNODAS model snow depth match an input raster") 
    parser.add_argument('-date', default=None, help='By default, date is extracted from input raster filename. Use this override or specify arbitrary timestamp (format: YYYYMMDD)')
    parser.add_argument('-datadir', default=os.getcwd(), help='Directory to store intermediate products (default: %(default)s)')
    parser.add_argument('fn', type=str, help='Raster filename to match (e.g., YYYYMMDD_raster.tif)')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    #Define top-level directory containing raster 
    topdir = os.getcwd()

    #This directory will store SNODAS products
    #Use centralized directory, default is $HOME/data/
    #datadir = iolib.get_datadir()
    datadir = args.datadir
    if not os.path.exists(datadir):
        os.makedirs(datadir)

    fn = args.fn
    ds = gdal.Open(fn)
    print(fn)

    #Extract timestamp from input filename
    dt = timelib.fn_getdatetime(fn)
    #If date is specified, extract timestamp 
    if args.date is not None:
        dt = timelib.fn_getdatetime(args.date)

    out_fn_base = os.path.splitext(fn)[0]

    snodas_min_dt = datetime(2003,9,30)
    if dt < snodas_min_dt:
        sys.exit("Timestamp is earlier than valid SNODAS model range")

    #snow depth values are mm, convert to meters
    snodas_outdir = os.path.join(datadir, 'snodas')
    if not os.path.exists(snodas_outdir):
        os.makedirs(snodas_outdir)
    snodas_ds_full = get_snodas(dt, snodas_outdir)
    snodas_ds = warplib.memwarp_multi([snodas_ds_full,], res='source', extent=ds, t_srs=ds, r='cubicspline')[0]
    snodas_depth = iolib.ds_getma(snodas_ds)/1000.

    if snodas_depth.count() > 0:
        #Write out at original resolution 
        out_fn = out_fn_base +'_snodas_depth.tif'
        print("Writing out %s" % out_fn)
        iolib.writeGTiff(snodas_depth, out_fn, src_ds=snodas_ds)

        #Warp to match input raster
        #Note: use cubicspline here to avoid artifacts with negative values
        ds_out = warplib.memwarp_multi([snodas_ds,], res=ds, extent=ds, t_srs=ds, r='cubicspline')[0]
        
        #Write out warped version
        snodas_depth = iolib.ds_getma(ds_out)/1000.
        out_fn = out_fn_base +'_snodas_depth_warp.tif'
        print("Writing out %s" % out_fn)
        iolib.writeGTiff(snodas_depth, out_fn, src_ds=ds_out)
    else:
        print("SNODAS grid for input location and timestamp is empty!")

if __name__ == "__main__":
    main()
