#!/bin/bash

today=`date +%Y%m`
prevMthStart=`date -d "$(date +%Y-%m-01) -1 day" +%Y-%m-01`
prevMthEnd=`date -d "$(date +%Y-%m-01) -1 day" +%Y-%m-%d`

ISAS_OA_DIR=/home/oo22/oo/co05/co0520/co052004/CORA-GLOBAL-04.1/OA/field/
ISAS_NRT_DIR=/home/coriolis_exp/spool/co04/co0401/ISAS_6_2/NRTOAGL01/ISAS_RESU/field/
WOA_DIR=/home/coriolis_dev/val/dat/co03/climatologie/WOA13/
CATDS_OPER_DIR=/home/catds_diss/data/cpdc/exp/current/FTP/OS/GRIDDED/L3OS/OPER/
CATDS_RE04_DIR=/home/catds_diss/data/cpdc/exp/current/FTP/OS/GRIDDED/L3OS/RE04/

echo `date` - Generating CATDS maps...
echo $prevMthStart
echo $prevMthEnd

# generate maps
cd /home/catds_sa/cpdc/vrec/install/current/v5.3/maps_catds/scripts/
./carte_comp.py -s $prevMthStart -e $prevMthEnd -o $ISAS_OA_DIR -n $ISAS_NRT_DIR -w $WOA_DIR -c $CATDS_OPER_DIR -r $CATDS_RE04_DIR

# copy maps to web server
cp *.png data/

# mkdir
mv_dir=png_$today
mkdir $mv_dir
mv *.png $mv_dir

echo `date` - End of generating CATDS maps.

