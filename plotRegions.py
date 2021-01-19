#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 11:09:39 2020

@author: sparky
"""

import shapefile
import numpy as np
import matplotlib as mpl
#mpl.use('AGG')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.path as mpltPath
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
import loadDF_winds as loadDF
import datetime
import pygrib
import os,sys
import pdb


# plots an image of all the subareas in each regional

rootfldr = os.path.abspath(os.path.dirname(sys.argv[0]))+'/'
shpfldr = rootfldr+'data/shapefiles/'
regionlist = ['NOTA','COTA','DUN','FLD','SN_LKS','SLD','CLU','MRB','NSN','WLD','BLR','CNY_HI','CPL','CH','NLD','ALD','BOP','COP','ROT','WKO','WTO','TPO','TKI','HBY','GSB','MNU','WGI','THP','TMN','WRP','KAP','WGN'] # Areas to do
nwpfn = '/data/models/ecmwf-prod/det/'
mdlrun = '2021011412'
NWPfn2 = rootfldr+'NWPfolder/NWPlatlongs.npy' # presaved lats and longs for domain
maplims = np.array([165.,186.,-33.,-48.]) # Total domain of NWP to load from
plt.close('all')
# Load lats and longs - SOME winds are differnet resolution
lats,lons = loadDF.loadLatLongs(False,NWPfn2,mdlrun,nwpfn,maplims) # Assume lats and longs are same through model run - res doesn't chane
# If res does change, just get lats and longs returned in loadNWP, and move spatial bit inside main time loop
# load in DFs

for region in regionlist:
    print('Processing '+region)
    
    shpfn = shpfldr + region+ 'subregions.shp'

# Load region, and clip data

    pts = np.array(zip(list(lons.flatten()), list(lats.flatten()))).astype('float')
    
    inArea = np.zeros(lats.shape)
    sf = shapefile.Reader(shpfn)
    
    pts = np.array(zip(list(lons.flatten()), list(lats.flatten()))).astype('float')
    sf = shapefile.Reader(shpfn)
    
    foundAll = False
    c = -1
    for shape in sf.shapeRecords():
        c +=1
        if shape.record[0] == 'All':
            foundAll = True
            if len(shape.shape.parts)==1:                           # Single part polygon
                i_start = shape.shape.parts[0]
                i_end = len(shape.shape.points)
                polygon = shape.shape.points[i_start:i_end]
                path = mpltPath.Path(polygon)
                inside = path.contains_points(pts)
                inArea = np.reshape(inside,lats.shape)
            else:                                                   # Deal with multiple parts to polygon one by one
                inSubTemp = np.zeros(lats.shape)
                for p in range(len(shape.shape.parts)):
                    i_start = shape.shape.parts[p]
                    if p<(len(shape.shape.parts)-1):
                        i_end = shape.shape.parts[p+1]
                    else:
                        i_end=len(shape.shape.points)
                    polygon = shape.shape.points[i_start:i_end]
                    path = mpltPath.Path(polygon)
                    inside = path.contains_points(pts)
                    inSubTemp += np.reshape(inside,lats.shape)
                inArea = inSubTemp>0
    
    if foundAll == False:
        print('No All region - check shapefile')
        sys.exit()
        
    # Load subregions
    nSDF = int(0) # number of subregions
    for shape in sf.shapeRecords():
        nSDF+=1
            
    SDF = np.zeros([lats.shape[0],lats.shape[1],nSDF]).astype(bool) # Dimensions lat, lon, subregion
    SDF_text = np.zeros(nSDF).astype('U256')
    
    c = -1    
    for shape in sf.shapeRecords():
        c+=1
        SDF_text[c] = shape.record[0]+' + '+shape.record[1]
        if shape.record[0] == 'All':
            wholeArea = c
        inverse = shape.record[1]
        if len(shape.shape.parts)==1:                           # Single part polygon
            i_start = shape.shape.parts[0]
            i_end = len(shape.shape.points)
            polygon = shape.shape.points[i_start:i_end]
            path = mpltPath.Path(polygon)
            inside = path.contains_points(pts)
            SDF[:,:,c] = np.reshape(inside,lats.shape)
        else:                                                   # Deal with multiple parts to polygon one by one
            inSubTemp = np.zeros(lats.shape)
            for p in range(len(shape.shape.parts)):
                i_start = shape.shape.parts[p]
                if p<(len(shape.shape.parts)-1):
                    i_end = shape.shape.parts[p+1]
                else:
                    i_end=len(shape.shape.points)
                polygon = shape.shape.points[i_start:i_end]
                path = mpltPath.Path(polygon)
                inside = path.contains_points(pts)
                inSubTemp += np.reshape(inside,lats.shape)
            SDF[:,:,c] = inSubTemp>0
        print('Subarea '+str(c)+' '+ SDF_text[c]+' '+str(sum(sum(SDF[:,:,c])))+' points')
    maplims2 = np.array([np.min(lons[inArea])-0.25,np.max(lons[inArea])+0.25,np.max(lats[inArea])+0.25,np.min(lats[inArea]-0.25)]) # Domain of region          
    fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
    nx = np.ceil(np.sqrt(nSDF)).astype(int)
    ny = np.ceil(nSDF/(1.0*nx)).astype(int)
    gs = mpl.gridspec.GridSpec(nx, ny)
    #gs.update(wspace=0.1, hspace=0.1, left=0.01, right=0.99, bottom=0.02, top=0.98) 
    
    for c in range(nSDF):
        ax=plt.subplot(gs[c])
        m = Basemap(llcrnrlon=maplims2[0],llcrnrlat=maplims2[3],urcrnrlon=maplims2[1],urcrnrlat=maplims2[2],projection='mill', resolution = 'h')
    #    step = np.diff(lons)[0,0]/2.0
    #    clats = lats+step
    #    clons = lons-step

    # Draw overall region
        shape = sf.shapeRecords()[wholeArea]
        i_start = shape.shape.parts[0]
        if len(shape.shape.parts)==1:                           # Single part polygon
            patches=[]        
            i_start = shape.shape.parts[0]
            i_end = len(shape.shape.points)
            rx = [i[0] for i in shape.shape.points[i_start:i_end]]
            ry = [i[1] for i in shape.shape.points[i_start:i_end]]
            x, y = m(rx, ry) 
            m.plot(x,y,'k-')
            poly = np.array(zip(list(x), list(y)))
            patches.append(mpatches.Polygon(poly))
            ax.add_collection(mpl.collections.PatchCollection(patches, facecolor= 'lightblue', edgecolor = 'k', linewidths =1.5))    
        # Draw regions
    
        shape = sf.shapeRecords()[c]
        i_start = shape.shape.parts[0]
        if len(shape.shape.parts)==1:                           # Single part polygon
            patches=[]        
            i_start = shape.shape.parts[0]
            i_end = len(shape.shape.points)
            rx = [i[0] for i in shape.shape.points[i_start:i_end]]
            ry = [i[1] for i in shape.shape.points[i_start:i_end]]
            x, y = m(rx, ry) 
            m.plot(x,y,'k-')
            poly = np.array(zip(list(x), list(y)))
            patches.append(mpatches.Polygon(poly))
            ax.add_collection(mpl.collections.PatchCollection(patches, facecolor= 'green', edgecolor = 'k', linewidths =1.5))
        else:                                                   # Deal with multiple parts to polygon one by one        
            patches=[]    
            for p in range(len(shape.shape.parts)):
                i_start = shape.shape.parts[p]
                if p<(len(shape.shape.parts)-1):
                    i_end = shape.shape.parts[p+1]
                else:
                    i_end=len(shape.shape.points)
                rx = [i[0] for i in shape.shape.points[i_start:i_end]]
                ry = [i[1] for i in shape.shape.points[i_start:i_end]]
                x, y = m(rx, ry)
                m.plot(x,y,'k-')
                poly = np.array(zip(list(x), list(y)))
                patches.append(mpatches.Polygon(poly))
                ax.add_collection(mpl.collections.PatchCollection(patches, facecolor= 'blue', edgecolor = 'k', linewidths =1.5))
    
        # Plot NWP points
        x, y = m(lons, lats)        
        m.plot(x,y,'k.',markersize=1)
        x, y = m(lons[SDF[:,:,c]], lats[SDF[:,:,c]])        
        m.plot(x,y,'go',markersize=5) 
        ax.set_title(SDF_text[c])
        m.drawcoastlines(color='k')
    fig.suptitle(region)
    imgfn = rootfldr+region+'.png'
    plt.savefig(imgfn)   
    