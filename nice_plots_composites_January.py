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
name_vortices_list = ['KOOBOR']#,'CAN2017_B2']#'CAN2017_A']#'CAN2017_A']
periods = [1, 2,1]
name_vortex = 'CAN2017_A'#'KOOBOR'
fig = plt.figure(figsize=(15.85,  1*5.))
First=True
bbox = dict(boxstyle='round4',fc='w')
ind_plot=['(a)','(b)','(c)', '(d)','(e)','(f)']

doincr = False
name='LPV'#Tdot'
def symbol(name):
    if (name == 'LPV'):
        return '$\\Pi$'
    else:
        return name

units = 'PVU'#' (K/day)'# ' (10$^{-6}$ kg/kg)'#' (PVU)'#'(10$^{-6}$ kg/kg)'#' (PVU)'
var = 'LPV'#'PV''LPV'#'O3'
doanom=False# True
def get_anomalie(aa, var='PV', doanom=True):
     field = aa[var].mean('time')
     if (doincr):
         field = (field-aa[var+'F'].mean('time'))/0.5
         if (var == 'T'):
             field = field + aa['PHR']*86400. # add heating rate to the increment
     elif (doanom):
         field = field-aa[var+'zon'].mean('time')
     if (('PV'  in var) or ('O3' in var)):
         field = 1.e6*field
     elif (var == 'PHR'):
         field = 86400.*field
     return field

exponent_LPV = 4.5#4.#4.5#9./2.
theta_ref = 420.
iki = 0
for ii, name_vortex in enumerate(name_vortices_list):
 nbase=131 +3*ii
 aa= xr.open_dataset('./interp_E5_inc_'+name_vortex+'_zz.nc')#interp_pressure_E5.nc') 
 if (name_vortex == '2ndVortex'):
     exponent_LPV = 4.75
 elif (name_vortex != 'KOOBOR'):
     exponent_LPV = 4.


 aa['PTzon'] = aa['Tzon']*(aa['P']/1.e5)**(-287./1004.)
 # introduce Lait PV to remove the vertical scaling
 aa['LPV'] = aa['PV']*(aa['PT']/theta_ref)**(-exponent_LPV)
 aa['LPVzon'] = aa['PVzon']*(aa['PTzon']/420.)**(-exponent_LPV)


 period=periods[ii]
 name_vortex_use = name_vortex.replace('_',' ') 

 if (name_vortex == 'KOOBOR'):
     period=periods[iki]
     iki=iki+1
     name_vortex_use = 'AUS2020 main'
     if (period==1):
         date1 = np.datetime64(datetime.datetime(2020,1,7))
         date2 = np.datetime64(datetime.datetime(2020,1,19))
     elif (period==2):
         date1 = np.datetime64(datetime.datetime(2020,2,2))
         date2 = np.datetime64(datetime.datetime(2020,2,28))
 elif ('CAN' in name_vortex):
         date1 = np.datetime64(datetime.datetime(2017,8,26))
         date2 = np.datetime64(datetime.datetime(2017,10,8))
 elif (name_vortex == '4thVortex'):
         date1 = np.datetime64(datetime.datetime(2020,1,15))
         date2 = np.datetime64(datetime.datetime(2020,2,15))
 elif (name_vortex == '2ndVortex'):
         date1 = np.datetime64(datetime.datetime(2020,1,18))
         date2 = np.datetime64(datetime.datetime(2020,1,28))


 
 
 
 aa = aa.sel(time=slice(date1,date2))
 date_str = pd.to_datetime(aa.time.values[0]).strftime('%Y/%m/%d') +' to '+pd.to_datetime(aa.time.values[-1]).strftime('%Y/%m/%d')
 ffac = 1.
 cmap = 'RdBu'
 if ('CAN' in name_vortex):
     ffac = 0.5
     if ('PV' in var):
         cmap = 'RdBu_r'
 if ('PHR' in var):
    cmap = 'RdBu_r'


 dmax=ffac*2.e3
 
 
 dz = 6.
 aa = aa.sel(xx=slice(-dmax,dmax)).sel(yy=slice(-dmax,dmax)).sel(zz=slice(-dz,dz))
 
 direc='./'
 
 field = get_anomalie(aa, var=var, doanom = doanom)
 VO = (aa['VO']).mean('time')
 
 seuils = np.array([3.e-5,5.e-5])
 if ('CAN' in name_vortex):
     seuils = ffac*np.sort(-seuils)

 levels = 20
 if (var == 'PHR'):
     levels = np.linspace(-0.8,0.8,17)

 ax2 = plt.subplot(nbase +1)
 coupe_horizontale = field.sel(yy=0., method='nearest').rename(var)
 coupe_hor_vo = VO.sel(yy=0., method='nearest').rename('VO')
 cf = ax2.contourf(coupe_horizontale.xx.values, coupe_horizontale.zz.values, coupe_horizontale.values, levels = levels, cmap = cmap, extend='both')
 ax2.contour(coupe_horizontale.xx.values, coupe_horizontale.zz.values,coupe_hor_vo.values,colors='g',levels = seuils)
 ax2.set_xlabel('x (km)')
 ax2.set_ylabel('z (km)')


 ax1 = plt.subplot(nbase)
 coupe_horizontale = field.sel(zz=0., method='nearest').rename(var)
 coupe_hor_vo = VO.sel(zz=0., method='nearest').rename('VO')
 cf = ax1.contourf(coupe_horizontale.xx.values, coupe_horizontale.yy.values, coupe_horizontale.values, levels = cf.levels, cmap = cmap, extend='both')
 ax1.contour(coupe_horizontale.xx.values, coupe_horizontale.yy.values,coupe_hor_vo.values,colors='g',levels = seuils)
 ax1.set_xlabel('x (km)')
 ax1.set_ylabel('y (km)')

 ffc = ffac
 xrepere =  coupe_horizontale.xx.values[0]+ffc*(coupe_horizontale.xx.values[16]-coupe_horizontale.xx.values[0])
 yrepere = coupe_horizontale.yy.values[-1] - ffc*(coupe_horizontale.yy.values[-1]-coupe_horizontale.yy.values[-4])
 
 
 ax1.text(xrepere,yrepere,name_vortex_use+' \n'+date_str, ha="center",va="center",size="13",bbox=bbox)

     

 ax3 = plt.subplot(nbase +2)
 coupe_horizontale = field.sel(xx=0., method='nearest').rename(var)
 coupe_hor_vo = VO.sel(xx=0., method='nearest').rename('VO')
 cf = ax3.contourf(coupe_horizontale.yy.values, coupe_horizontale.zz.values, coupe_horizontale.values, levels = cf.levels, cmap = cmap, extend='both')
 ax3.contour(coupe_horizontale.yy.values, coupe_horizontale.zz.values,coupe_hor_vo.values,colors='g',levels = seuils)
 ax3.set_xlabel('y (km)')
 ax3.set_ylabel('z (km)')
 
 levs = list(np.linspace(-5.,5.,11))
 levs= ffac*np.array(levs)
 levs = levs[np.abs(levs)>=1.] 
 field = get_anomalie(aa, var='T')
 
 coupe_horizontale = field.sel(yy=0., method='nearest').rename(var)
 CS = ax2.contour(coupe_horizontale.xx.values, coupe_horizontale.zz.values, coupe_horizontale.values, levels = levs, colors='k')
 ax2.clabel(CS, inline=1, fontsize=10, fmt = '%1.1f')
 coupe_horizontale = field.sel(zz=0., method='nearest').rename(var)
 CS = ax1.contour(coupe_horizontale.xx.values, coupe_horizontale.yy.values, coupe_horizontale.values, levels = levs, colors='k')
 ax1.clabel(CS, inline=1, fontsize=10, fmt = '%1.1f')
 
 coupe_horizontale = field.sel(xx=0., method='nearest').rename(var)
 CS = ax3.contour(coupe_horizontale.yy.values, coupe_horizontale.zz.values, coupe_horizontale.values, levels = levs,colors='k')
 ax3.clabel(CS, inline=1, fontsize=10, fmt = '%1.1f')

 pos = np.array(matplotlib.artist.getp(ax3, 'position'))
 cbar_ax = fig.add_axes([0.92, pos[0,1], 0.02, pos[1,1]-pos[0,1]])
 cbar = plt.colorbar(cf, cax=cbar_ax, orientation = 'vertical')
 cbar.set_label(symbol(var)+ ' ('+units+')')

 ypos=1.06
 if (First):
     First= False
     ax1.set_title('horizontal section', y=ypos)
     ax2.set_title('zonal section', y=ypos)
     ax3.set_title('meridional section', y=ypos)
 for ij, ax in enumerate([ax1, ax2, ax3]): 
     ik = 3*ii + ij
     ax.text(-0.1, 1.03, ind_plot[ik], transform=ax.transAxes, size=14, weight='bold')
#name='Tdot'
plt.savefig('./figures_papier/nice_composites_'+var+'_'+str(period)+'_January.pdf')#,dpi=300)
plt.show()



'''
plt.suptitle(var+' increment '+name_vortex)
plt.figure()
field.mean('time').sel(xx=0., method='nearest').rename(var).drop('xx').contour.contourf(levels = 20)
(aa['VO']).mean('time').sel(xx=0., method='nearest').rename('VO').drop('xx').contour.contour(colors='k',levels = seuils)
plt.title(var+' increment '+name_vortex)
plt.savefig(direc+'composite_vertmeridian_increment'+name_vortex+'_'+var+'.png')

plt.figure()
field.mean('time').sel(yy=0., method='nearest').rename(var).drop('yy').contour.contourf(levels = 20)
(aa['VO']).mean('time').sel(yy=0., method='nearest').rename('VO').drop('yy').contour.contour(colors='k',levels = seuils)
plt.title(var+' increment '+name_vortex)
plt.savefig(direc+'composite_vertzonal_increment'+name_vortex+'_'+var+'.png')

plt.show()
'''
