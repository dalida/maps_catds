#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: tabstop=8 expandtab softtabstop=4 shiftwidth=4 textwidth=79
import os
import glob
import tarfile
import math
import datetime
import numpy as np
import scipy
import scipy.interpolate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import netCDF4
import logging
from optparse import OptionParser

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


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

# parsing options
    parser = OptionParser()
    parser.add_option("-s", "--start", action = "store", type = "string", dest = "start", default = '2010-01-01')
    parser.add_option("-e", "--end", action = "store", type = "string", dest = "end", default = datetime.datetime.now().strftime('%Y-%m-01'))
    (options, args) = parser.parse_args()

    try:
        start = datetime.datetime.strptime(options.start,'%Y-%m-%d')
    except AttributeError:
        print 'Error start'
    try:
        end = datetime.datetime.strptime(options.end,'%Y-%m-%d')
    except AttributeError:
        print 'Error end'
      

    date = start
    while date<end:
        logging.info('Date : {}'.format(date))

        isas_exists = False
        isas_dir = '/net/argos/data/ipso/DATA/MyOCEAN/INSITU_GLO_TS_OA_REP_OBSERVATIONS_013_002_b/CORIOLIS-GLOBAL-CORA04.1-OBS_FULL_TIME_SERIE/field/{}/'.format(date.year)
        isas_file = 'OA_CORA4.1_{}15_fld_PSAL.nc' .format(date.strftime('%Y%m'))
        if not(os.path.exists(isas_dir + isas_file)):
            isas_dir = '/net/argos/data/ipso/DATA/MyOCEAN/INSITU_GLO_TS_OA_NRT_OBSERVATIONS_013_002_a/CORIOLIS-GLOBAL-NRTOA-OBS/ISAS_RESU/field/%04d/' % ( date.year )
            isas_file = 'OA_NRTOAGL01_{}15_fld_PSAL.nc'.format(date.strftime('%Y%m'))
        if os.path.exists(isas_dir + isas_file):
            isas_exists = True
            isas_nc = netCDF4.Dataset (isas_dir + isas_file)
            logging.info(isas_file + ' .... loaded !')
            isas_lat = isas_nc.variables['latitude'][:]
            isas_lon = isas_nc.variables['longitude'][:]
            isas_sss = isas_nc.variables['PSAL'][0, 0, :, :] #2nd field : prof 0=>0m 1=>3m 2=>5m
            isas_pcv = isas_nc.variables['PCTVAR'][0, 0, :, :]
            isas_nc.close()
            isas_sss = np.ma.masked_array(isas_sss, mask=isas_pcv > 80, fill_value=np.nan)
            isas_mlon, isas_mlat = np.meshgrid(isas_lon, isas_lat)

        if isas_exists:
            isas_title = 'ISAS [{}]'.format(date.strftime('%Y-%m'))
            title = isas_title
            isas_fig = 'ISAS_{}'.format(date.strftime('%Y%m'))
            fig_file = isas_fig + '.png'

            fig = plt.figure()
            ax=fig.add_subplot(111)
            ax.set_title(title)
            plt.figtext(1,0,isas_file,ha='right',va='bottom',fontsize=6)
            map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
            map.drawcoastlines(linewidth=0.01,antialiased=False)
            map.drawmapboundary(fill_color='white',linewidth=0.01)
            map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
            map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
            map.fillcontinents(color='grey')
            N=12
            vmin=32
            vmax=38
            ticks = np.linspace(vmin,vmax,N+1)
            isas_lonout,z =map.shiftdata(isas_lon, isas_sss, lon_0=-50)
            lon,lat = np.meshgrid(isas_lonout,isas_lat)
            x,y=map(lon,lat)
            cmap = cm.get_cmap('jet',N)
            cmap.set_bad('1.0')
            cmap.set_under((0.0, 0.0, 0.25, 1.0))
            cmap.set_over((0.25, 0.0, 0.0, 1.0))
            pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
            cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
            cb.set_label('SSS [PSS]')
            plt.savefig(fig_file)
            logging.info(fig_file +' .... created !' )
            plt.close()

        woa_exists = False
        woa_dir = '/net/argos/data/ipso/DATA/WOA/2013/'
        woa_file = 'woa13_A5B2_s{:02d}_04.nc'.format(date.month)
        if os.path.exists(woa_dir + woa_file):
            woa_exists = True
            try:
                woa_nc = netCDF4.Dataset (woa_dir + woa_file)
                woa_lat = woa_nc.variables['lat'][:]
                woa_lon = woa_nc.variables['lon'][:]
                woa_sss = woa_nc.variables['s_an'][0,0,:,:]
                woa_nc.close()
                logging.info(woa_file + ' .... loaded !')
            except RuntimeError:
                continue

            index = np.argsort(((woa_lon + 180)%360) -180)
            woa_lon = ((woa_lon[index]+180.)%360.)-180.
            woa_sss = woa_sss[:,index]
            woa_sss = np.ma.masked_array(woa_sss, mask=(woa_sss<0.), fill_value=np.nan)
            woa_mlon, woa_mlat = np.meshgrid(woa_lon, woa_lat)
        else:
            logging.error('{} not exists !'.format(woa_dir + woa_file))

        if woa_exists:
            woa_title = 'World Ocean Atlas $_{2013}$ '+'[{}]'.format(date.strftime('%m'))
            title = woa_title
            woa_fig = 'WOA2013_{}'.format(date.strftime('%m'))
            fig_file = woa_fig + '.png'

            fig = plt.figure()
            ax=fig.add_subplot(111)
            ax.set_title(title)
            plt.figtext(1,0,woa_file,ha='right',va='bottom',fontsize=6)
            map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
            map.drawcoastlines(linewidth=0.01,antialiased=False)
            map.drawmapboundary(fill_color='white',linewidth=0.01)
            map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
            map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None)
            map.fillcontinents(color='grey')
            N=12
            vmin=32
            vmax=38
            ticks = np.linspace(vmin,vmax,N+1)
            woa_lonout,z =map.shiftdata(woa_lon, woa_sss, lon_0=-50)
            lon,lat = np.meshgrid(woa_lonout,woa_lat)
            x,y=map(lon,lat)
            cmap = cm.get_cmap('jet',N)
            cmap.set_bad('1.0')
            cmap.set_under((0.0, 0.0, 0.25, 1.0))
            cmap.set_over((0.25, 0.0, 0.0, 1.0))
            pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
            cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
            cb.set_label('SSS [PSS]')
            plt.savefig(fig_file)
            logging.info(fig_file +' .... created !' )
            plt.close()

        if woa_exists and isas_exists:
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
                mask=isas2woa_mask, 
                fill_value=np.nan,
                )

            title = isas_title + ' - ' + woa_title
            fig_file = isas_fig + '-' + woa_fig + '.png'

            fig = plt.figure()
            ax=fig.add_subplot(111)
            ax.set_title(title)
            map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
            map.drawcoastlines(linewidth=0.01,antialiased=False)
            map.drawmapboundary(fill_color='white',linewidth=0.01)
            map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
            map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
            map.fillcontinents(color='grey')
            N=10
            vmin=-0.5
            vmax=+0.5
            ticks = np.linspace(vmin,vmax,N+1)
            woa_lonout,z =map.shiftdata(woa_lon, isas2woa_sss - woa_sss, lon_0=-50)
            lon,lat = np.meshgrid(woa_lonout,woa_lat)
            x,y=map(lon,lat)
            cmap = cm.get_cmap('RdBu_r',N)
            cmap.set_bad('1.0')
            cmap.set_under((0.0, 0.0, 0.2, 1.0))
            cmap.set_over((0.2, 0.0, 0.0, 1.0))
            pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
            cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
            cb.set_label('SSS [PSS]')
            plt.savefig(fig_file)
            logging.info(fig_file +' .... created !' )
            plt.close()

        for orbit in ('A','D','AD',):
            smosLOCEAN_exists = False
            smosLOCEAN_dirs = [
                    '/net/argos/data/ipso/DATA/SMOS/L3_LOCEANv2013_REPR_L2OS_v622/MONTH/',
                    '/net/argos/data/ipso/DATA/SMOS/L3_LOCEANv2013/MONTH/',
                    ]
            conds = [
                    'REPR', 
                    'OCOR', 
                    'OPER', 
                    ]
            for smosLOCEAN_dir in smosLOCEAN_dirs:
                for cond in conds :
                    smosLOCEAN_file = 'SMOS_L3_LOCEAN_%s_%s_%04d%02d01-%04d%02d01_0.25deg_100.0km.nc' % (cond, orbit, date.year, date.month, date.year + date.month/12, (date.month+1-1)%12+1 )
                    if os.path.exists(smosLOCEAN_dir + smosLOCEAN_file):
                        smosLOCEAN_exists = True
                        break
            if smosLOCEAN_exists:
                logging.info(smosLOCEAN_file + ' .... loaded !')
                smosLOCEAN_nc = netCDF4.Dataset(smosLOCEAN_dir + smosLOCEAN_file)
                smosLOCEAN_lat = smosLOCEAN_nc.variables['lat'][:]
                smosLOCEAN_lon = smosLOCEAN_nc.variables['lon'][:]
                smosLOCEAN_sss = smosLOCEAN_nc.variables['SSS'][0, :, :]
                smosLOCEAN_n = smosLOCEAN_nc.variables['nSSS'][0, :, :]
                smosLOCEAN_nc.close()
                smosLOCEAN_lon = np.fmod(smosLOCEAN_lon+180.+3600.0, 360.)-180.
                index = np.argsort(smosLOCEAN_lon)
                smosLOCEAN_lon = smosLOCEAN_lon[index].flatten()
                smosLOCEAN_lat = smosLOCEAN_lat.flatten()
                smosLOCEAN_sss = smosLOCEAN_sss[:, index]
                smosLOCEAN_n = smosLOCEAN_n[:, index]
                mask = np.ma.mask_or(
                    np.ma.make_mask(smosLOCEAN_sss > 50.),
                    np.ma.make_mask(smosLOCEAN_n < 30.)
                    )
                smosLOCEAN_sss = np.ma.masked_array(
                    smosLOCEAN_sss,
                    mask=mask,
                    fill_value=np.nan,
                    )
                smosLOCEAN_mlon, smosLOCEAN_mlat = np.meshgrid(smosLOCEAN_lon, smosLOCEAN_lat)

                smosLOCEAN_title = 'N SMOS LOCEAN-CEC [{},{}]'.format(orbit,date.strftime('%Y-%m'))
                title = smosLOCEAN_title
                smosLOCEAN_fig = 'N_SMOS_LOCEAN_{}_{}'.format(orbit,date.strftime('%Y%m'))
                fig_file = smosLOCEAN_fig + '.png'

                fig = plt.figure()
                ax=fig.add_subplot(111)
                ax.set_title(title)
                plt.figtext(1,0,smosLOCEAN_file,ha='right',va='bottom',fontsize=6)
                map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
                map.drawcoastlines(linewidth=0.01,antialiased=False)
                map.drawmapboundary(fill_color='white',linewidth=0.01)
                map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
                map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
                map.fillcontinents(color='grey')
                N=15
                vmin=0
                vmax=300
                ticks = np.linspace(vmin,vmax,N+1)
                smosLOCEAN_lonout,z =map.shiftdata(smosLOCEAN_lon, smosLOCEAN_n, lon_0=-50)
                lon,lat = np.meshgrid(smosLOCEAN_lonout,smosLOCEAN_lat)
                x,y=map(lon,lat)
                cmap = cm.get_cmap('jet',N)
                cmap.set_bad('1.0')
                cmap.set_under((0.0, 0.0, 0.2, 1.0))
                cmap.set_over((0.2, 0.0, 0.0, 1.0))
                pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
                cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'max',ticks=ticks)
                cb.set_label('SSS [PSS]')
                plt.savefig(fig_file)
                logging.info(fig_file +' .... created !' )
                plt.close()


                smosLOCEAN_title = 'SMOS LOCEAN-CEC[{},{}]'.format(orbit,date.strftime('%Y-%m'))
                title = smosLOCEAN_title
                smosLOCEAN_fig = 'SMOS_LOCEAN_{}_{}'.format(orbit,date.strftime('%Y%m'))
                fig_file = smosLOCEAN_fig + '.png'

                fig = plt.figure()
                ax=fig.add_subplot(111)
                ax.set_title(title)
                plt.figtext(1,0,smosLOCEAN_file,ha='right',va='bottom',fontsize=6)
                map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
                map.drawcoastlines(linewidth=0.01,antialiased=False)
                map.drawmapboundary(fill_color='white',linewidth=0.01)
                map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
                map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
                map.fillcontinents(color='grey')
                N=12
                vmin=32
                vmax=38
                ticks = np.linspace(vmin,vmax,N+1)
                smosLOCEAN_lonout,z =map.shiftdata(smosLOCEAN_lon, smosLOCEAN_sss, lon_0=-50)
                lon,lat = np.meshgrid(smosLOCEAN_lonout,smosLOCEAN_lat)
                x,y=map(lon,lat)
                cmap = cm.get_cmap('jet',N)
                cmap.set_bad('1.0')
                cmap.set_under((0.0, 0.0, 0.2, 1.0))
                cmap.set_over((0.2, 0.0, 0.0, 1.0))
                pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
                cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
                cb.set_label('SSS [PSS]')
                plt.savefig(fig_file)
                logging.info(fig_file +' .... created !' )
                plt.close()
                
                date0 = date
                date1 = (date0+datetime.timedelta(days=45)).replace(day=1)-datetime.timedelta(seconds=1)
                smosCATDS_exists = False
                smosCATDS_dir = '/net/argos/data/ipso/DATA/SMOS/CATDS-CPDC/OPER/MIR_CSF3B{}/{:04d}/{:02d}/'.format(
                        '_' if orbit=='AD' else orbit,
                        date0.year,
                        date0.month,)
# SM_OPER_MIR_CSF3B__20150801T000000_20150831T235959_300_002_7.tgz
                smosCATDS_file = 'SM_OPER_MIR_CSF3B{}_{}_{}_*_*_7.tgz'.format(
                        '_' if orbit=='AD' else orbit,
                        date0.strftime('%Y%m%dT%H%M%S'),
                        date1.strftime('%Y%m%dT%H%M%S'),
                        )
                smosCATDS_files = glob.glob(smosCATDS_dir + smosCATDS_file)
                if len(smosCATDS_files)>0:
                    smosCATDS_file = os.path.basename(smosCATDS_files[-1])
                if os.path.exists(smosCATDS_dir + smosCATDS_file):
                    smosCATDS_exists = True
                if smosCATDS_exists:
                    tf = tarfile.open(smosCATDS_dir + smosCATDS_file)
                    for tm in tf.getmembers():

                        if tm.name[-3:]=='.nc':
                            tf.extract(tm,path='/tmp')
                            break

                    logging.info('/tmp/'+tm.name + ' .... loaded !')
                    smosCATDS_nc = netCDF4.Dataset ('/tmp/' + tm.name)
                    smosCATDS_lat = smosCATDS_nc.variables['lat'][:]
                    smosCATDS_lon = smosCATDS_nc.variables['lon'][:]
                    smosCATDS_sss = smosCATDS_nc.variables['Mean_Sea_Surface_Salinity'][:,:]
                    smosCATDS_n = smosCATDS_nc.variables['N_Used_Meas'][:,:]
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
                        mask=mask,
                        fill_value=np.nan,
                        )
                    smosCATDS_n = np.ma.masked_array(
                        smosCATDS_n,
                        mask=np.ma.make_mask(smosCATDS_n==0),
                        fill_value=0,
                        )
                    smosCATDS_mlon, smosCATDS_mlat = np.meshgrid(smosCATDS_lon, smosCATDS_lat)


                    smosCATDS_title = 'N SMOS CATDS-CPDC [{},{}]'.format(orbit,date.strftime('%Y-%m'))
                    title = smosCATDS_title
                    smosCATDS_fig = 'N_SMOS_CATDS_{}_{}'.format(orbit,date.strftime('%Y%m'))
                    fig_file = smosCATDS_fig + '.png'

                    fig = plt.figure()
                    ax=fig.add_subplot(111)
                    ax.set_title(title)
                    plt.figtext(1,0,smosCATDS_file,ha='right',va='bottom',fontsize=6)
                    map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
                    map.drawcoastlines(linewidth=0.01,antialiased=False)
                    map.drawmapboundary(fill_color='white',linewidth=0.01)
                    map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
                    map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
                    map.fillcontinents(color='grey')
                    N=15
                    vmin=0
                    vmax=75
                    ticks = np.linspace(vmin,vmax,N+1)
                    smosCATDS_lonout,z =map.shiftdata(smosCATDS_lon, smosCATDS_n, lon_0=-50)
                    lon,lat = np.meshgrid(smosCATDS_lonout,smosCATDS_lat)
                    x,y=map(lon,lat)
                    cmap = cm.get_cmap('jet',N)
                    cmap.set_bad('1.0')
                    cmap.set_under((0.0, 0.0, 0.2, 1.0))
                    cmap.set_over((0.2, 0.0, 0.0, 1.0))
                    pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
                    cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'max',ticks=ticks)
                    cb.set_label('SSS [PSS]')
                    plt.savefig(fig_file)
                    logging.info(fig_file +' .... created !' )
                    plt.close()

                    smosCATDS_title = 'SMOS CATDS-CPDC [{},{}]'.format(orbit,date.strftime('%Y-%m'))
                    title = smosCATDS_title
                    smosCATDS_fig = 'SMOS_CATDS_{}_{}'.format(orbit,date.strftime('%Y%m'))
                    fig_file = smosCATDS_fig + '.png'

                    fig = plt.figure()
                    ax=fig.add_subplot(111)
                    ax.set_title(title)
                    plt.figtext(1,0,smosCATDS_file,ha='right',va='bottom',fontsize=6)
                    map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
                    map.drawcoastlines(linewidth=0.01,antialiased=False)
                    map.drawmapboundary(fill_color='white',linewidth=0.01)
                    map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
                    map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
                    map.fillcontinents(color='grey')
                    N=12
                    vmin=32
                    vmax=38
                    ticks = np.linspace(vmin,vmax,N+1)
                    smosCATDS_lonout,z =map.shiftdata(smosCATDS_lon, smosCATDS_sss, lon_0=-50)
                    lon,lat = np.meshgrid(smosCATDS_lonout,smosCATDS_lat)
                    x,y=map(lon,lat)
                    cmap = cm.get_cmap('jet',N)
                    cmap.set_bad('1.0')
                    cmap.set_under((0.0, 0.0, 0.2, 1.0))
                    cmap.set_over((0.2, 0.0, 0.0, 1.0))
                    pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
                    cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
                    cb.set_label('SSS [PSS]')
                    plt.savefig(fig_file)
                    logging.info(fig_file +' .... created !' )
                    plt.close()
            
                if smosCATDS_exists:
                    smosCATDS2smosLOCEAN_sss = scipy.interpolate.griddata(
                        (smosCATDS_mlon.flatten(), smosCATDS_mlat.flatten()),
                        smosCATDS_sss.flatten(),
                        (smosLOCEAN_mlon, smosLOCEAN_mlat),
                        method='linear',
                        )
                    smosCATDS2smosLOCEAN_mask = scipy.interpolate.griddata(
                        (smosCATDS_mlon.flatten(),smosCATDS_mlat.flatten()),
                        smosCATDS_sss.mask.flatten(),
                        (smosLOCEAN_mlon, smosLOCEAN_mlat), 
                        method='linear',
                        )
                    smosCATDS2smosLOCEAN_sss = np.ma.masked_array(
                        smosCATDS2smosLOCEAN_sss,
                        mask=smosCATDS2smosLOCEAN_mask, 
                        fill_value=np.nan,
                        )

                if isas_exists:
                    isas2smosLOCEAN_sss = scipy.interpolate.griddata(
                        (isas_mlon.flatten(), isas_mlat.flatten()),
                        isas_sss.flatten(),
                        (smosLOCEAN_mlon, smosLOCEAN_mlat),
                        method='linear',
                        )
                    isas2smosLOCEAN_mask = scipy.interpolate.griddata(
                        (isas_mlon.flatten(),isas_mlat.flatten()),
                        isas_sss.mask.flatten(),
                        (smosLOCEAN_mlon, smosLOCEAN_mlat), 
                        method='linear',
                        )
                    isas2smosLOCEAN_sss = np.ma.masked_array(
                        isas2smosLOCEAN_sss,
                        mask=isas2smosLOCEAN_mask, 
                        fill_value=np.nan,
                        )
                if woa_exists:
                    woa2smosLOCEAN_sss = scipy.interpolate.griddata(
                        (woa_mlon.flatten(), woa_mlat.flatten()),
                        woa_sss.flatten(),
                        (smosLOCEAN_mlon, smosLOCEAN_mlat),
                        method='linear',
                        )
                    woa2smosLOCEAN_mask = scipy.interpolate.griddata(
                        (woa_mlon.flatten(),woa_mlat.flatten()),
                        woa_sss.mask.flatten(),
                        (smosLOCEAN_mlon, smosLOCEAN_mlat), 
                        method='linear',
                        )
                    woa2smosLOCEAN_sss = np.ma.masked_array(
                        woa2smosLOCEAN_sss,
                        mask=woa2smosLOCEAN_mask, 
                        fill_value=np.nan,
                        )
        

                if smosCATDS_exists:
                    title = smosCATDS_title + ' - ' + smosLOCEAN_title
                    fig_file = smosCATDS_fig + '-' + smosLOCEAN_fig + '.png'

                    fig = plt.figure()
                    ax=fig.add_subplot(111)
                    ax.set_title(title)
                    map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
                    map.drawcoastlines(linewidth=0.01,antialiased=False)
                    map.drawmapboundary(fill_color='white',linewidth=0.01)
                    map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
                    map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
                    map.fillcontinents(color='grey')
                    N=10
                    vmin=-0.5
                    vmax=+0.5
                    ticks = np.linspace(vmin,vmax,N+1)
                    smosLOCEAN_lonout,z =map.shiftdata(smosLOCEAN_lon, smosCATDS2smosLOCEAN_sss -  smosLOCEAN_sss, lon_0=-50)
                    lon,lat = np.meshgrid(smosLOCEAN_lonout,smosLOCEAN_lat)
                    x,y=map(lon,lat)
                    cmap = cm.get_cmap('RdBu_r',N)
                    cmap.set_bad('1.0')
                    cmap.set_under((0.0, 0.0, 0.2, 1.0))
                    cmap.set_over((0.2, 0.0, 0.0, 1.0))
                    pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
                    cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
                    cb.set_label('SSS [PSS]')
                    plt.savefig(fig_file)
                    logging.info(fig_file +' .... created !' )
                    plt.close()
                    

                if isas_exists:
                    title = smosLOCEAN_title + ' - ' + isas_title
                    fig_file = smosLOCEAN_fig + '-' + isas_fig + '.png'

                    fig = plt.figure()
                    ax=fig.add_subplot(111)
                    ax.set_title(title)
                    map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
                    map.drawcoastlines(linewidth=0.01,antialiased=False)
                    map.drawmapboundary(fill_color='white',linewidth=0.01)
                    map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
                    map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
                    map.fillcontinents(color='grey')
                    N=10
                    vmin=-0.5
                    vmax=+0.5
                    ticks = np.linspace(vmin,vmax,N+1)
                    smosLOCEAN_lonout,z =map.shiftdata(smosLOCEAN_lon, smosLOCEAN_sss - isas2smosLOCEAN_sss, lon_0=-50)
                    lon,lat = np.meshgrid(smosLOCEAN_lonout,smosLOCEAN_lat)
                    x,y=map(lon,lat)
                    cmap = cm.get_cmap('RdBu_r',N)
                    cmap.set_bad('1.0')
                    cmap.set_under((0.0, 0.0, 0.2, 1.0))
                    cmap.set_over((0.2, 0.0, 0.0, 1.0))
                    pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
                    cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
                    cb.set_label('SSS [PSS]')
                    plt.savefig(fig_file)
                    logging.info(fig_file +' .... created !' )
                    plt.close()
                    
                    title = smosCATDS_title + ' - ' + isas_title
                    fig_file = smosCATDS_fig + '-' + isas_fig + '.png'

                    fig = plt.figure()
                    ax=fig.add_subplot(111)
                    ax.set_title(title)
                    map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
                    map.drawcoastlines(linewidth=0.01,antialiased=False)
                    map.drawmapboundary(fill_color='white',linewidth=0.01)
                    map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
                    map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
                    map.fillcontinents(color='grey')
                    N=10
                    vmin=-0.5
                    vmax=+0.5
                    ticks = np.linspace(vmin,vmax,N+1)
                    smosLOCEAN_lonout,z =map.shiftdata(smosLOCEAN_lon, smosCATDS2smosLOCEAN_sss - isas2smosLOCEAN_sss, lon_0=-50)
                    lon,lat = np.meshgrid(smosLOCEAN_lonout,smosLOCEAN_lat)
                    x,y=map(lon,lat)
                    cmap = cm.get_cmap('RdBu_r',N)
                    cmap.set_bad('1.0')
                    cmap.set_under((0.0, 0.0, 0.2, 1.0))
                    cmap.set_over((0.2, 0.0, 0.0, 1.0))
                    pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
                    cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
                    cb.set_label('SSS [PSS]')
                    plt.savefig(fig_file)
                    logging.info(fig_file +' .... created !' )
                    plt.close()
                    

                if woa_exists:
                    title = smosLOCEAN_title + ' - ' + woa_title
                    fig_file = smosLOCEAN_fig + '-' + woa_fig + '.png'

                    fig = plt.figure()
                    ax=fig.add_subplot(111)
                    ax.set_title(title)
                    map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
                    map.drawcoastlines(linewidth=0.01,antialiased=False)
                    map.drawmapboundary(fill_color='white',linewidth=0.01)
                    map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
                    map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
                    map.fillcontinents(color='grey')
                    N=10
                    vmin=-0.5
                    vmax=+0.5
                    ticks = np.linspace(vmin,vmax,N+1)
                    smosLOCEAN_lonout,z =map.shiftdata(smosLOCEAN_lon, smosLOCEAN_sss - woa2smosLOCEAN_sss, lon_0=-50)
                    lon,lat = np.meshgrid(smosLOCEAN_lonout,smosLOCEAN_lat)
                    x,y=map(lon,lat)
                    cmap = cm.get_cmap('RdBu_r',N)
                    cmap.set_bad('1.0')
                    cmap.set_under((0.0, 0.0, 0.2, 1.0))
                    cmap.set_over((0.2, 0.0, 0.0, 1.0))
                    pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
                    cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
                    cb.set_label('SSS [PSS]')
                    plt.savefig(fig_file)
                    logging.info(fig_file +' .... created !' )
                    plt.close()
                    
                    title = smosCATDS_title + ' - ' + woa_title
                    fig_file = smosCATDS_fig + '-' + woa_fig + '.png'

                    fig = plt.figure()
                    ax=fig.add_subplot(111)
                    ax.set_title(title)
                    map = Basemap(projection='moll',resolution='l',lon_0=-50,ellps='WGS84',anchor='S')
                    map.drawcoastlines(linewidth=0.01,antialiased=False)
                    map.drawmapboundary(fill_color='white',linewidth=0.01)
                    map.drawmeridians(np.arange(-180,181,60), labels=[0,0,0,0],linewidth=0.01,labelstyle=None)
                    map.drawparallels(np.arange(-90,91,30), labels=[1,0,0,0],linewidth=0.01,labelstyle=None) 
                    map.fillcontinents(color='grey')
                    N=10
                    vmin=-0.5
                    vmax=+0.5
                    ticks = np.linspace(vmin,vmax,N+1)
                    smosLOCEAN_lonout,z =map.shiftdata(smosLOCEAN_lon, smosCATDS2smosLOCEAN_sss - woa2smosLOCEAN_sss, lon_0=-50)
                    lon,lat = np.meshgrid(smosLOCEAN_lonout,smosLOCEAN_lat)
                    x,y=map(lon,lat)
                    cmap = cm.get_cmap('RdBu_r',N)
                    cmap.set_bad('1.0')
                    cmap.set_under((0.0, 0.0, 0.2, 1.0))
                    cmap.set_over((0.2, 0.0, 0.0, 1.0))
                    pc=map.pcolormesh(x,y,z,vmin=vmin,vmax=vmax,cmap=cmap)
                    cb=plt.colorbar(pc,shrink=0.8,orientation='horizontal',fraction=0.04, extend = 'both',ticks=ticks)
                    cb.set_label('SSS [PSS]')
                    plt.savefig(fig_file)
                    logging.info(fig_file +' .... created !' )
                    plt.close()

        date = datetime.datetime( date.year + date.month/12, (date.month+1-1)%12+1, 1,0,0,0,0,None)
