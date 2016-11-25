#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: tabstop=8 expandtab softtabstop=4 shiftwidth=4 textwidth=79

'''
    carte_comp.py
    Generates CATDS, ISAS, WOA maps and the comparisons between them.

    2016-10-16 LM Initial

    2016-11-15 LM WOA and CATDS directories as parameters

'''

import os
import glob
import tarfile
import math
import datetime
import numpy as np
import scipy
import scipy.interpolate
import netCDF4
#import ncdump
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
from optparse import OptionParser

# configure plot
matplotlib.rc('figure' ,autolayout='True')
matplotlib.rc('lines', antialiased='False')
matplotlib.rc('patch', antialiased='False')
matplotlib.rc('font', family='sans-serif',size=8, style='normal',weight='normal')
matplotlib.rc('text', hinting='auto',usetex='False')
matplotlib.rc('mathtext', default='sf')
matplotlib.rc('axes', linewidth=1,titlesize=8,labelsize=6)
matplotlib.rc('xtick', labelsize=6)
matplotlib.rc('ytick', labelsize=6)
matplotlib.rc('figure', dpi=100, figsize=(12,6))
matplotlib.rc('savefig', dpi=100)

def createMap(title_in, file_in, fig_file_in, N, vmin, vmax, lon_in,
             lat_in, sss_in, colors, label='SSS [PSS]'):
    """
    Creates a map with given input parameters.
    """

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_title(title_in)
    plt.figtext(1, 0, file_in, ha='right', va='bottom', fontsize=6)

    map = Basemap(projection='moll', resolution='l', lon_0=-50, ellps='WGS84', anchor='S')
    map.drawcoastlines(linewidth=0.01, antialiased=False)
    map.drawmapboundary(fill_color='white', linewidth=0.01)
    map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0], linewidth=0.01, labelstyle=None)
    map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0], linewidth=0.01, labelstyle=None) 
    map.fillcontinents(color='grey')

    ticks = np.linspace(vmin, vmax, N+1)
    
    lonout, z = map.shiftdata(lon_in, sss_in, lon_0=-50)
    lon, lat  = np.meshgrid(lonout, lat_in)
    x, y      = map(lon, lat)

    cmap = cm.get_cmap(colors, N)
    cmap.set_bad('1.0')
    cmap.set_under((0.0, 0.0, 0.25, 1.0))
    cmap.set_over((0.25, 0.0, 0.0, 1.0))

    pc = map.pcolormesh(x, y, z, vmin=vmin, vmax=vmax, cmap=cmap)
    cb = plt.colorbar(pc, shrink=0.8, orientation='horizontal', fraction=0.04, extend ='both', ticks=ticks)
    cb.set_label(label)
    plt.savefig(fig_file_in)
    logging.debug(fig_file_in +' .... created!' )
    plt.close()

    return None



if __name__ == '__main__':

    # config logging
    logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

    # parsing options
    parser = OptionParser()
    parser.add_option("-s", "--start", action = "store", type = "string", dest = "start", default = '2010-01-01')
    parser.add_option("-e", "--end", action = "store", type = "string", dest = "end", default = datetime.datetime.now().strftime('%Y-%m-01'))
    parser.add_option("-c", "--catds_dir", action = "store", type = "string", dest = "catds_dir")
    parser.add_option("-w", "--woa_dir", action = "store", type = "string",
dest = "woa_dir")
    (options, args) = parser.parse_args()

    # define dates to process
    try:
        start = datetime.datetime.strptime(options.start,'%Y-%m-%d')
    except AttributeError:
        print 'Error with option start'
        log.info('Error with option start')
    try:
        end = datetime.datetime.strptime(options.end,'%Y-%m-%d')
    except AttributeError:
        print 'Error with option end'
        log.info('Error with option end')
    try:
        catds_dir = options.catds_dir 
    except AttributeError:
        print 'Error with option catds_dir'
        log.info('Error with option catds_dir')
    try:
        woa_dir = options.woa_dir 
    except AttributeError:
        print 'Error with option woa_dir'
        log.info('Error with option woa_dir')
      

    date = start
    while date < end:
        logging.info('start - {}'.format(start) + ' - end - {}'.format(end))

        # Generate ISAS maps
        isas_exists = False

        isas_dir  = '/home/oo22/oo/co05/co0520/co052004/CORA-GLOBAL-04.1/OA/field/%04d/' % (date.year )
        isas_file = 'OA_CORA4.1_{}15_fld_PSAL.nc'.format(date.strftime('%Y%m'))

        logging.debug('Reprocessed isas path {} ...'.format(isas_dir + isas_file))

        if not(os.path.exists(isas_dir + isas_file)):
            logging.debug('isas_file does not exist...')
            isas_dir = '/home/coriolis_exp/spool/co04/co0401/ISAS_6_2/NRTOAGL01/ISAS_RESU/field/%04d/' % ( date.year )
            isas_file = 'OA_NRTOAGL01_{}15_fld_PSAL.nc'.format(date.strftime('%Y%m'))

        if os.path.exists(isas_dir + isas_file):
            logging.debug('isas_file exists...')
            isas_exists = True


        # Generate ISAS maps
        if isas_exists:
            logging.info('  Generating ISAS maps...')

            isas_nc = netCDF4.Dataset (isas_dir + isas_file)
            logging.debug(isas_file + ' .... netcdf file loaded!')

            # DEBUG
            #ncdump.ncdump(isas_nc)

            # get variables from netcdf files
            isas_lat = isas_nc.variables['latitude'][:]
            isas_lon = isas_nc.variables['longitude'][:]
            isas_sss = isas_nc.variables['PSAL'][0, 0, :, :] #2nd field : prof 0=>0m 1=>3m 2=>5m
            isas_pcv = isas_nc.variables['PCTVAR'][0, 0, :, :]
            isas_nc.close()
            isas_sss = np.ma.masked_array(isas_sss, mask=isas_pcv > 80, fill_value=np.nan)
            isas_mlon, isas_mlat = np.meshgrid(isas_lon, isas_lat)

            isas_title = 'ISAS [{}]'.format(date.strftime('%Y-%m'))
            title = isas_title
            isas_fig = 'ISAS_{}'.format(date.strftime('%Y%m'))
            fig_file = isas_fig + '.png'

            createMap(title, isas_file, fig_file, 12, 32, 38, isas_lon,
                     isas_lat, isas_sss, 'jet')

        # Generate WOA maps
        #woa_dir = '/home/coriolis_dev/val/dat/co03/climatologie/WOA13/'
        woa_exists = False
        woa_file = 's{:02d}_04.nc'.format(date.month)

        if os.path.exists(woa_dir + woa_file):
            woa_exists = True

        # Generate WOA maps
        if woa_exists:
            logging.info('  Generating WOA maps...')
            try:
                woa_nc = netCDF4.Dataset (woa_dir + woa_file)
                woa_lat = woa_nc.variables['lat'][:]
                woa_lon = woa_nc.variables['lon'][:]
                woa_sss = woa_nc.variables['s_an'][0,0,:,:]
                #logging.debug(woa_file + ' .... loaded !')
            except RuntimeError:
                continue

            index = np.argsort(((woa_lon + 180)%360) -180)
            woa_lon = ((woa_lon[index]+180.)%360.)-180.
            woa_sss = woa_sss[:,index]
            woa_sss = np.ma.masked_array(woa_sss, mask=(woa_sss<0.), fill_value=np.nan)
            woa_mlon, woa_mlat = np.meshgrid(woa_lon, woa_lat)

            woa_title = 'World Ocean Atlas $_{2013}$ '+'[{}]'.format(date.strftime('%m'))
            woa_title = 'World Ocean Atlas 2013 '+'[{}]'.format(date.strftime('%m'))
            title = woa_title
            woa_fig = 'WOA2013_{}'.format(date.strftime('%m'))
            fig_file = woa_fig + '.png'

            createMap(title, woa_file, fig_file, 11, 32, 38, woa_lon,
                       woa_lat, woa_sss, 'jet')

        else:
            logging.error('woa dir and file {} does not exist !'.format(woa_dir + woa_file))


        # Generate ISAS and WOA comparison map 
        if isas_exists and woa_exists:
            logging.info('  Generating WOA and ISAS comparison map...')

            isas2woa_sss = scipy.interpolate.griddata(
                (isas_mlon.flatten(), isas_mlat.flatten()),
                isas_sss.flatten(),
                (woa_mlon, woa_mlat),
                method='linear',
                )
            isas2woa_mask = scipy.interpolate.griddata(
                (isas_mlon.flatten(),isas_mlat.flatten()),
                isas_sss.mask.flatten(),
                (woa_mlon, woa_mlat), 
                method='linear',
                )
            isas2woa_sss = np.ma.masked_array(
                isas2woa_sss,
                mask = isas2woa_mask, 
                fill_value=np.nan,
                )

            isas_woa_sss = isas2woa_sss - woa_sss

            title = isas_title + ' - ' + woa_title
            fig_file = isas_fig + '-' + woa_fig + '.png'

            createMap(title, woa_file, fig_file, 10, -0.5, +0.5, woa_lon,
                       woa_lat, isas_woa_sss, 'RdBu_r')

        # flags for A-D
        resetA, resetD = 0, 0

        # Generate CATDS maps
        for orbit in ('A','D','AD',):

            date0 = date
            date1 = (date0+datetime.timedelta(days=45)).replace(day=1)-datetime.timedelta(seconds=1)
            smosCATDS_exists = False

            #smosCATDS_dir = '/home/catds_diss/data/cpdc/exp/current/FTP/OS/GRIDDED/L3OS/OPER/MIR_CSF3B{}/{:04d}/{:02d}/'.format('_'if orbit=='AD' else orbit,
            smosCATDS_dir = catds_dir + 'MIR_CSF3B{}/{:04d}/{:02d}/'.format('_'if orbit=='AD' else orbit,
                date0.year,
                date0.month,) 
            # currently only 2015 and 2016
            smosCATDS_file = 'SM_OPER_MIR_CSF3B{}_{}_{}_*_*_7.tgz'.format('_' if orbit=='AD' else orbit,
                date0.strftime('%Y%m%dT%H%M%S'),
                date1.strftime('%Y%m%dT%H%M%S'),)
            logging.debug('CATDS dir ' + smosCATDS_dir)
            logging.debug('Processing file ' + smosCATDS_file)

            smosCATDS_files = glob.glob(smosCATDS_dir + smosCATDS_file)
            if len(smosCATDS_files)>0:
                smosCATDS_file = os.path.basename(smosCATDS_files[-1])
            if os.path.exists(smosCATDS_dir + smosCATDS_file):
                smosCATDS_exists = True

            if smosCATDS_exists:
                logging.info('  Generating CATDS maps')
                logging.debug('smosCATDS exists...')
                tf = tarfile.open(smosCATDS_dir + smosCATDS_file)
                for tm in tf.getmembers():
                    if tm.name[-3:]=='.nc':
                        tf.extract(tm,path='/tmp')
                        break

                logging.debug('/tmp/'+tm.name + ' .... loaded !')
                smosCATDS_nc = netCDF4.Dataset ('/tmp/' + tm.name)
                smosCATDS_lat = smosCATDS_nc.variables['lat'][:]
                smosCATDS_lon = smosCATDS_nc.variables['lon'][:]
                smosCATDS_sss = smosCATDS_nc.variables['Mean_Sea_Surface_Salinity'][:,:]
                smosCATDS_n = smosCATDS_nc.variables['N_Used_Meas'][:,:]
                smosCATDS_err = smosCATDS_nc.variables['Sss_Standard_Deviation'][:,:]
                smosCATDS_std = smosCATDS_nc.variables['Sss_Rms_Mean'][:,:]
                smosCATDS_nc.close()
                smosCATDS_lon = np.fmod(smosCATDS_lon+180.+3600.0, 360.)-180.
                index = np.argsort(smosCATDS_lon)
                smosCATDS_lon = smosCATDS_lon[index].flatten()
                smosCATDS_lat = smosCATDS_lat.flatten()
                smosCATDS_sss = smosCATDS_sss[:, index]
                smosCATDS_n = smosCATDS_n[:, index]
                mask = np.ma.mask_or(
                    np.ma.make_mask(smosCATDS_sss > 50.),
                    np.ma.make_mask(smosCATDS_n < 8.)
                )
                smosCATDS_sss = np.ma.masked_array(
                    smosCATDS_sss,
                    mask = mask,
                    fill_value=np.nan,
                )
                smosCATDS_n = np.ma.masked_array(
                    smosCATDS_n,
                    mask = np.ma.make_mask(smosCATDS_n==0),
                    fill_value=0,
                )
                smosCATDS_err = np.ma.masked_array(
                    smosCATDS_err,
                    mask=mask,
                    fill_value=np.nan,
                )
                smosCATDS_std = np.ma.masked_array(
                    smosCATDS_std,
                    mask=mask,
                    fill_value=np.nan,
                )

                smosCATDS_mlon, smosCATDS_mlat = np.meshgrid(smosCATDS_lon, smosCATDS_lat)


                # N SMOS CATDS
                smosCATDS_title = 'N SMOS CATDS_CPDC [{},{}]'.format(orbit,date.strftime('%Y-%m'))
                title = smosCATDS_title
                smosCATDS_fig = 'N_SMOS_CATDS_{}_{}'.format(orbit,date.strftime('%Y%m'))
                fig_file = smosCATDS_fig + '.png'

                createMap(title, smosCATDS_file, fig_file, 15, 0, 75, smosCATDS_lon,
                            smosCATDS_lat, smosCATDS_n, 'jet', 'Number of Used Measures')

                # SMOS CATDS
                smosCATDS_title = 'SMOS CATDS-CPDC [{},{}]'.format(orbit,date.strftime('%Y-%m'))
                title = smosCATDS_title
                smosCATDS_fig = 'SMOS_CATDS_{}_{}'.format(orbit,date.strftime('%Y%m'))
                fig_file = smosCATDS_fig + '.png'

                createMap(title, smosCATDS_file, fig_file, 12, 32, 38, smosCATDS_lon,
                            smosCATDS_lat, smosCATDS_sss, 'jet')

                # STD/sqrt(N) (Sss_Standard_Deviation in NetCDF file)
                smosCATDSerr_title = 'Std/sqrt(N)' + ' SMOS CATDS_CPDC [{},{}]'.format(orbit,date.strftime('%Y-%m'))
                title = smosCATDSerr_title
                smosCATDSerr_fig = 'ERR_SMOS_CATDS_{}_{}'.format(orbit,date.strftime('%Y%m'))
                fig_file = smosCATDSerr_fig + '.png'

                createMap(title, smosCATDS_file, fig_file, 12, 0, +0.6, smosCATDS_lon,
                            smosCATDS_lat, smosCATDS_err, 'jet', 'Std(SSS)/sqrt(N) [PSS]')

                # STD (Sss_Rms_Mean in NetCDF file)
                smosCATDSstd_title = 'Std SMOS CATDS_CPDC [{},{}]'.format(orbit,date.strftime('%Y-%m'))
                title = smosCATDSstd_title
                smosCATDSstd_fig = 'STD_SMOS_CATDS_{}_{}'.format(orbit,date.strftime('%Y%m'))
                fig_file = smosCATDSstd_fig + '.png'

                createMap(title, smosCATDS_file, fig_file, 10, 0, +2, smosCATDS_lon,
                            smosCATDS_lat, smosCATDS_std, 'jet', 'Std(SSS) [PSS]')

	        # generate A-D map if we've got both A and D
                # A-D
	        if orbit=='A':
	            resetA = 1
		    smosa_file = smosCATDS_file
		    smosa_sss = smosCATDS_sss
		    smosa_path = smosCATDS_exists
	        elif orbit=='D':
	            resetD = 1
		    smosd_file = smosCATDS_file
		    smosd_sss = smosCATDS_sss
		    smosd_path = smosCATDS_exists
		    smosd_lon = smosCATDS_lon
	            smosd_lat = smosCATDS_lat
	        else:
		    smosa_file, smosd_file = None, None
		    smosa_sss, smosd_sss = None, None
		    smosd_lon, smosd_lat = None, None
		    resetA, resetD = 0, 0
                
	        if resetA==1 and resetD==1 and smosa_path and smosd_path:

                    #smosad_lonout,z =map.shiftdata(smosd_lon, smosa_sss - smosd_sss, lon_0=-50)
                    #lon,lat = np.meshgrid(smosad_lonout,smosd_lat)
		    sss_diff = smosa_sss - smosd_sss
		    smosad_file = smosa_file + ' + ' + smosd_file

                    smosad_title = 'SMOS CATDS_CPDC [A-D,{}]'.format(date.strftime('%Y-%m'))
                    title = smosad_title
                    smosad_fig = 'SMOS_CATDS_AmD_{}'.format(date.strftime('%Y%m'))
                    fig_file = smosad_fig + '.png'

	            createMap(title, smosad_file, fig_file, 10, -0.5, +0.5, smosd_lon, smosd_lat,
		        sss_diff, 'RdBu_r', 'delta SSS [PSS]')

                # CATDS - ISAS comparison map
                if isas_exists:
                    isas2catds_sss = scipy.interpolate.griddata(
                        (isas_mlon.flatten(), isas_mlat.flatten()),
                        isas_sss.flatten(),
                        (smosCATDS_mlon, smosCATDS_mlat),
                        method='linear',
                        )
                    isas2catds_mask = scipy.interpolate.griddata(
                        (isas_mlon.flatten(),isas_mlat.flatten()),
                        isas_sss.mask.flatten(),
                        (smosCATDS_mlon, smosCATDS_mlat), 
                        method='linear',
                        )
                    isas2catds_sss = np.ma.masked_array(
                        isas2catds_sss,
                        mask = isas2catds_mask, 
                        fill_value=np.nan,
                    )
                
                    isas_catds_sss = smosCATDS_sss - isas2catds_sss

                    title = smosCATDS_title + ' - ' + isas_title
                    fig_file = smosCATDS_fig + '-' + isas_fig + '.png'

                    createMap(title, smosCATDS_file, fig_file, 10, -0.5, +0.5, smosCATDS_lon,
                               smosCATDS_lat, isas_catds_sss, 'RdBu_r')

                # CATDS - WOA comparison map
                if woa_exists:
                    woa2catds_sss = scipy.interpolate.griddata(
                        (woa_mlon.flatten(), woa_mlat.flatten()),
                        woa_sss.flatten(),
                        (smosCATDS_mlon, smosCATDS_mlat),
                        method='linear',
                        )
                    woa2catds_mask = scipy.interpolate.griddata(
                        (woa_mlon.flatten(),woa_mlat.flatten()),
                        woa_sss.mask.flatten(),
                        (smosCATDS_mlon, smosCATDS_mlat), 
                        method='linear',
                        )
                    woa2catds_sss = np.ma.masked_array(
                        woa2catds_sss,
                        mask = woa2catds_mask, 
                        fill_value=np.nan,
                    )
                
                    woa_catds_sss = smosCATDS_sss - woa2catds_sss

                    title = smosCATDS_title + ' - ' + woa_title
                    fig_file = smosCATDS_fig + '-' + woa_fig + '.png'

                    createMap(title, smosCATDS_file, fig_file, 10, -0.5, +0.5, smosCATDS_lon,
                               smosCATDS_lat, woa_catds_sss, 'RdBu_r')

        # Increment to next date to process
        date = datetime.datetime( date.year + date.month/12, (date.month+1-1)%12+1, 1,0,0,0,0,None)
