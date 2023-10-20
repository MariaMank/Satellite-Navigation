# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 17:48:22 2021

@author: Maciek
"""

#import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc, rcParams
import matplotlib
import numpy as np
import math as mat
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, 'D:/Maciek/PW/python')
sys.path.insert(1, 'E:/ZGIAG/SNS/satpos_almanac')
#from positioning import date2tow
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def latlon(XYZ):
    """
    function which converts ECEF satellite position position to latitude and longitude
    Based on rvtolatlong.m in Richard Rieber's orbital library on mathwork.com

    Note: that this is only suitable for groundtrack visualization, not rigorous
    calculations.
    """
    lati = []
    longti = []
    a = len(XYZ[0])
    for i  in range(a):
        xyz = [XYZ[0][i], XYZ[1][i], XYZ[2][i]]
        r_delta = np.linalg.norm(xyz[0:1]);

        sinA = xyz[1]/r_delta;
        cosA = xyz[0]/r_delta;

        Lon = mat.atan2(sinA,cosA)
        if Lon < - mat.pi:
            Lon = Lon + 2 * mat.pi;
        Lat = mat.asin(xyz[2]/np.linalg.norm(xyz))
        lati.append(mat.degrees(Lat))
        longti.append(mat.degrees(Lon))

    return lati, longti



# sats_data = {'1':np.array(([[82.3842,	-91.9837,	1.61529e+07,	3.34739e+06,	2.06534e+07],[84.1748,	-132.482,	1.67285e+07,	4.80881e+06,	1.9922e+07]]))}
# # sats_data = np.load('my_file.npy',allow_pickle='TRUE').item()
# time = np.array(([[86400, 0],[87000,	600]]))
# data = [2021,3,1,0,0,0]
# sats_data = {'1':np.array(([[82.3842,	-91.9837,	1.61529e+07,	3.34739e+06,	2.06534e+07],[84.1748,	-132.482,	1.67285e+07,	4.80881e+06,	1.9922e+07]]))}

# sats_to_plot = 'all'; MASKA = 10; visible= True
#%%
def groundtrack(sat_data, BLr, MASKA):
        # a = 6378137
    # def groundtrack(sats_data):
    # rc('text', usetex=True)
    # rc('font', family='serif');
    # rc('font',family='helvetica');
    rc('grid', color='gray', linewidth=0.5, linestyle='--')
    fontsize = 20
    rc('xtick', labelsize = fontsize)
    rc('ytick', labelsize = fontsize)
    rc('font', size = fontsize)
    fig =plt.figure(figsize=(14,7))
    plt.subplots_adjust(bottom= 0.1, 
                            top   = 0.9,
                            left  = 0.1, 
                            right = 0.9)
    ax = fig.add_subplot()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()
    ax.grid(True)
    plt.yticks(range(-90,100,30));
    plt.xticks([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180]);
    for sat in sats_data:
        data_one_sat = sats_data[sat][:,0:]
        bl_one_sat = np.zeros((0,2))
        for data_one_epoch in data_one_sat:
            bl_sat_one_epoch = np.array(([latlon(data_one_epoch[2:])]))
            bl_one_sat = np.append(bl_one_sat,bl_sat_one_epoch, axis=0)
        
        ax.plot(bl_one_sat[:,1],bl_one_sat[:,0],color='b')

    plt.show()

if __name__ == '__main__':
    sats_data = {'1':np.array(([[82.3842,	-91.9837,	1.61529e+07,	3.34739e+06,	2.06534e+07],[84.1748,	-132.482,	1.67285e+07,	4.80881e+06,	1.9922e+07]]))}
    # sats_data = np.load('my_file.npy',allow_pickle='TRUE').item()
    time = np.array(([[86400, 0],[87000,	600]]))
    data = [2021,3,1,0,0,0]
    groundtrack(sats_data, [52,21], 15) 