#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 13:54:18 2020

shows the composite

@author: aurelien
"""

import numpy as np
import gzip,pickle
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
name_vortices_list = ['KOOBOR','CAN2017_A']#,'CAN2017_B1']#'CAN2017_A']
periods = [2,1]
suffix = 'KOOBOR'
name_vortex = 'CAN2017_A'#'KOOBOR'
fig = plt.figure(figsize=(15.85,  2*5.))
First=True
bbox = dict(boxstyle='round4',fc='w')#, background="w")
ind_plot=['(a)','(b)','(c)','(d)','(e)','(f)']

doincr = True
units = 'PVU$\\, $ d$^{-1}$'
var = 'T'
if (var =='T'):
    units = 'K$\\, $ d$^{-1}$'

def get_anomalie(aa, var='PV', doincr = doincr, doanom=True):
     field = aa[var]
     if (doincr):
         field = (field-(aa[var+'F']))/0.5
         if (var == 'T'):
             field = field #+ aa['PHR']*86400. # add heating rate to the increment
     elif (doanom):
         field = (field-aa[var+'zon'])
     field = field.mean('time')
     if (var=='PV'):
         field = 1.e6*field
     elif (var == 'PHR'):
         field = 86400.*field
     return field

def slanted_axis(xcoord, alpha):
    return xcoord*np.cos(alpha), xcoord*np.sin(alpha)


for ii, name_vortex in enumerate(name_vortices_list):
 nbase=221 +2*ii
 aa= xr.open_dataset('./interp_E5_inc_'+name_vortex+'_zz.nc')#interp_pressure_E5.nc') 

 if (name_vortex == '2ndVortex'):
     exponent_LPV = 4.75
     theta_ref = 420.
 elif (name_vortex != 'KOOBOR'):
     exponent_LPV = 4.
     theta_ref = 420.

 lev_cols = 20
 if (var =='PV'):
     lev_cols = np.linspace(-6.,6., 17)
 elif (var=='T'):
     lev_cols = np.linspace(-.8,.8,17)
 period=periods[ii]
 alpha=0.
 name_vortex_use = name_vortex.replace('_',' ') 
 if (name_vortex == 'KOOBOR'):
     name_vortex_use = 'AUS2020 main'
     alpha = 0.#np.pi/4.-np.pi/2.
     if (period==1):
         date1 = np.datetime64(datetime.datetime(2020,1,7))
         date2 = np.datetime64(datetime.datetime(2020,1,20))
     elif (period==2):
         date1 = np.datetime64(datetime.datetime(2020,2,2))
         date2 = np.datetime64(datetime.datetime(2020,2,28))
 if ('CAN' in name_vortex):
         date1 = np.datetime64(datetime.datetime(2017,8,26))
         date2 = np.datetime64(datetime.datetime(2017,10,8))
         if (var == 'PV'):
             lev_cols = np.linspace(-2.,2., 17)
         elif (var=='T'):
             lev_cols = np.linspace(-0.3,0.3,13)
 
 aa = aa.sel(time=slice(date1,date2))
 date_str = pd.to_datetime(aa.time.values[0]).strftime('%Y/%m/%d') +' to '+pd.to_datetime(aa.time.values[-1]).strftime('%Y/%m/%d')

 ffac = 1.
 cmap = 'RdBu'
 if ('CAN' in name_vortex):
     ffac = 0.5
     if (var == 'PV'):
         cmap = 'RdBu_r'

 dz = 8.
 if (var=='T'):
     cmap = 'RdBu_r'
 dmax=ffac*2.e3
 dz = ffac*dz 




 aa = aa.sel(xx=slice(-dmax,dmax)).sel(yy=slice(-dmax,dmax)).sel(zz=slice(-dz,dz))
 
 direc='./'
 
 field = get_anomalie(aa, var=var)
 field_PHR = get_anomalie(aa, var='PHR', doincr=False, doanom=False)

 VO = (aa['VO']).mean('time')
 
 seuils = np.array([3.e-5,5.e-5])
 if ('CAN' in name_vortex):
     seuils = ffac*np.sort(-seuils)
 xx_a, yy_a = slanted_axis(field.xx.values, alpha)
 xx_b, yy_b = xx_a, yy_a#slanted_axis(field.xx.values, alpha+np.pi/2.)


 ax3 = plt.subplot(nbase+1 )
 coupe_horizontale = field.interp(xx=('xhor', xx_a), yy=('xhor', yy_a)).rename(var)
 coupe_hor_vo = VO.interp(xx=('xhor', xx_a), yy=('xhor', yy_a)).rename('VO')
 cf = ax3.contourf(coupe_horizontale.xx.values, coupe_horizontale.zz.values, coupe_horizontale.values, levels = lev_cols, cmap = cmap, extend='both')
 ax3.contour(coupe_horizontale.xx.values, coupe_horizontale.zz.values,coupe_hor_vo.values,colors='g',levels = seuils)
 xlab = 'x$_a$ (km)' if (np.abs(alpha)>1.e-2) else 'x (km)'
 ax3.set_xlabel(xlab)
 ax3.set_ylabel('z (km)')
 if (First):
     xrepere = 1./ffac*coupe_horizontale.xx.values[11]
     yrepere = 1./ffac*coupe_horizontale.zz.values[-3]
 else:
     xrepere = ffac*xrepere
     yrepere = ffac*yrepere

     
 ax2 = plt.subplot(nbase)


 coupe_horizontale_PHR = field_PHR.interp(xx=('xhor', xx_b), yy=('xhor', yy_b)).rename(var)
 coupe_hor_vo = VO.interp(xx=('xhor', xx_b), yy=('xhor', yy_b)).rename('VO')

 '''
 coupe_horizontale = field.sel(xx=0., method='nearest').rename(var)
 coupe_hor_vo = VO.sel(xx=0., method='nearest').rename('VO')
 '''
 cf = ax2.contourf(coupe_horizontale_PHR.xx.values, coupe_horizontale_PHR.zz.values, coupe_horizontale_PHR.values, levels = lev_cols, cmap = cmap, extend='both')
 ax2.contour(coupe_horizontale.xx.values, coupe_horizontale.zz.values,coupe_hor_vo.values,colors='g',levels = seuils)
 xlab = 'x$_b$ (km)' if (np.abs(alpha)>1.e-2) else 'x (km)'
 ax2.set_xlabel(xlab)
 ax2.set_ylabel('z (km)')



 levs = list(np.linspace(-5.,5.,11))
 levs= ffac*np.array(levs)
 levs = levs[np.abs(levs)>=1.]
 field = get_anomalie(aa, var='T', doincr = False)

 if ((var == 'T') or var=='PV'):
     coupe_horizontale = field.sel(yy=0., method='nearest').rename(var)
     CS = ax3.contour(coupe_horizontale.xx.values, coupe_horizontale.zz.values, coupe_horizontale.values, levels = levs, colors='k')
     ax3.clabel(CS, inline=1, fontsize=10, fmt = '%1.1f')

     CS = ax2.contour(coupe_horizontale.xx.values, coupe_horizontale.zz.values, coupe_horizontale.values, levels = levs,colors='k')
     ax2.clabel(CS, fontsize=10, fmt = '%1.1f')#, inline=1)
 pos = np.array(matplotlib.artist.getp(ax2, 'position'))
 cbar_ax = fig.add_axes([0.92, pos[0,1], 0.02, pos[1,1]-pos[0,1]])
 cbar = plt.colorbar(cf, cax=cbar_ax, orientation = 'vertical')
 cbar.set_label(units)

 ax2.text(xrepere,yrepere,name_vortex_use+' \n'+date_str, ha="center",va="center",size="13",bbox=bbox, zorder=1000000)



 ypos=1.06
 if (First):
     First= False
     if (np.abs(alpha)>1.e-2):
         xtitle = 'x$_a$ section'
         ytitle = 'y$_a$ section'
     else:
         xtitle = 'heating rate'
         ytitle = 'T increment tendency'

     ax2.set_title(xtitle, y=ypos)
     ax3.set_title(ytitle, y=ypos)
 for ij, ax in enumerate([ax2, ax3]): 
     ik = 2*ii + ij
     ax.text(-0.1, 1.03, ind_plot[ik], transform=ax.transAxes, size=14, weight='bold')

plt.savefig('./figures_papier/nice_composites_heating_rate.pdf')#,dpi=300)
plt.show()



