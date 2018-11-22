#! /bin/bash

#Create mosaics and snow depth products after DEM generation and co-registration

#Test case for 2017 Grand Mesa WV3 products 

res=8
#res=2

outdir=dem_nomap_mos_${res}m

#Want to cluster DEMs acquired within a few days
#See timelib utilities

#For now, hardcode
#ts_list="20160925 20170127 20170201 20170226 2017031[78]"
ts_list="20160925 20170127 20170201 20170226"

parallel "dem_mosaic -o $outdir/{}_mos_${res}m.tif *{}*00/dem_nomap/*-DEM_${res}m.tif" ::: $ts_list
dem_mosaic -o $outdir/20170318_mos_${res}m.tif *20170317*00/dem_nomap/*-DEM_${res}m.tif *20170318*00/dem_nomap/*-DEM_${res}m.tif 

cd $outdir
#Compute all elevation difference combinations
all_dh.sh *_mos.tif

#Generate composites
dem_mosaic -o all_snowdepth_${res}m_20170127-20170318.tif 20160925_mos_${res}m*eul.tif
dem_mosaic -o 20160925_mos_${res}m_20170127-20170201_mos_dz_eul.tif 20160925_mos_${res}m_20170127_mos_dz_eul.tif 20160925_mos_${res}m_20170201_mos_dz_eul.tif

