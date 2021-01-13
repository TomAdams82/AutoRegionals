#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 22:08:46 2021

@author: tadams@met.co.nz
"""
import numpy as np
import matplotlib as mpl
#mpl.use('AGG')
import matplotlib.pyplot as plt
#from matplotlib import gridspec
#import matplotlib.path as mpltPath
#import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
import pdb

def plotTimeSeriesMaps(data_timestep,lats,lons,TP,sf,maplims,inArea,TDF,SDF_text,SDF,SSSn,txt):
    fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
    gs = mpl.gridspec.GridSpec(3, 3)
    #gs.update(wspace=0.1, hspace=0.1, left=0.01, right=0.99, bottom=0.02, top=0.98) 
    # Find record for All Area (not the same as wholeArea which relates to inverses too)
    c  = -1
    for shape in sf.shapeRecords():
        c +=1
        if shape.record[0] == 'All':
            wholeAreaRecord = c
        if shape.record[0] == SDF_text:
            subAreaRecord = c
        if shape.record[1] == SDF_text:
            subAreaRecord = c
            
    for t in range(24/data_timestep):
        ax=plt.subplot(gs[t]) 
        m = Basemap(llcrnrlon=maplims[0],llcrnrlat=maplims[3],urcrnrlon=maplims[1],urcrnrlat=maplims[2],projection='mill', resolution = 'h')
    #    step = np.diff(lons)[0,0]/2.0
    #    clats = lats+step
    #    clons = lons-step
        m.drawcoastlines(color='g')
    # Draw overall region
        shape = sf.shapeRecords()[wholeAreaRecord]
        i_start = shape.shape.parts[0]
        if len(shape.shape.parts)==1:                           # Single part polygon
            i_start = shape.shape.parts[0]
            i_end = len(shape.shape.points)
            rx = [i[0] for i in shape.shape.points[i_start:i_end]]
            ry = [i[1] for i in shape.shape.points[i_start:i_end]]
            x, y = m(rx, ry) 
            m.plot(x,y,'k-')
        else:                                                   # Deal with multiple parts to polygon one by one        
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
        
        # Plot precip
        step = np.diff(lons)[0,0]/2.0
        clats = lats+step
        clons = lons-step
        x, y = m(clons, clats)        
        #m.pcolor(x, y, tp, norm=colors.LogNorm(vmin=0, vmax=1), cmap='Blues')
        #    m.pcolor(x, y, tp, vmin=0, vmax=0.5, cmap='Blues')    
        m.pcolormesh(x, y, TP[:,:,t], vmin=0, vmax=10, cmap='Blues')#
        m.colorbar()
        #If has a spatial extent, draw regions
        if ((TDF[t]==0) & (SDF_text != 'All')):
            shape = sf.shapeRecords()[subAreaRecord]
            i_start = shape.shape.parts[0]
            if len(shape.shape.parts)==1:                           # Single part polygon
                i_start = shape.shape.parts[0]
                i_end = len(shape.shape.points)
                rx = [i[0] for i in shape.shape.points[i_start:i_end]]
                ry = [i[1] for i in shape.shape.points[i_start:i_end]]
                x, y = m(rx, ry) 
                m.plot(x,y,'k-')
            else:                                                   # Deal with multiple parts to polygon one by one        
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
        # Plot NWP points
        x, y = m(lons, lats)        
        m.plot(x,y,'k.',markersize=1)
        utc_labels = np.concatenate([(np.arange(15,24,3)),(np.arange(0,15,3))]).astype(str)
        ax.set_title(utc_labels[t]+'Z',fontsize=14)
        if ((TDF[t]==0) & (SDF_text != 'All')):
            subregion_YYY = SDF[:,:,SSSn]
            subregion_XXX = inArea * np.logical_not(subregion_YYY)
            x, y = m(lons[subregion_XXX], lats[subregion_XXX])        
            m.plot(x,y,'go',markersize=5) 
            x, y = m(lons[subregion_YYY], lats[subregion_YYY])        
            m.plot(x,y,'rx',markersize=5)
        else:
            subregion_XXX = inArea
            x, y = m(lons[subregion_XXX], lats[subregion_XXX])        
            if (TDF[t]==0):
                m.plot(x,y,'go',markersize=5) 
            else:
                m.plot(x,y,'rx',markersize=5)
    fig.suptitle(txt)

def plotTimeSeriesDFs(TP_bins,nwpDF_XXX,nwpDF_YYY,guessDF_XXX,guessDF_YYY,scoreXXX,scoreYYY,data_timestep,txt):  
    fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
    gs = mpl.gridspec.GridSpec(3, 3)
    nlabels = len(TP_bins)-1  
    precip_labels = np.zeros(nlabels+2).astype(str)
    precip_labels[0]=''
    precip_labels[-1]=''
    utc_labels = np.concatenate([(np.arange(15,24,3)),(np.arange(0,15,3))]).astype(str)
    for i in range(nlabels):
        precip_labels[i+1] = str(TP_bins[i])+'-'+str(TP_bins[i+1])
    x = np.arange(nlabels)
    for t in range(24/data_timestep):
        #score[i] = KL_divergence(nwpDF[:,t], guessDF[:,t], base=None) 
        
        ax=plt.subplot(gs[t])    
        ax.plot(x,nwpDF_XXX[:,t],'b-')
        ax.plot(x,guessDF_XXX[:,t],'b:')
        ax.plot(x,nwpDF_YYY[:,t],'r-')
        ax.plot(x,guessDF_YYY[:,t],'r:')
        ax.set_title(utc_labels[t]+'Z, KL-div score '+str(scoreXXX[t])+' XXX, '+str(scoreYYY[t])+' YYY',fontsize=14)
        #ax.set_yticks(x)
        ax.set_xlim([-1,6])
        ax.set_xticklabels(precip_labels)
        ax.set_xlabel('Precip amount / mm')
        ax.set_ylim([0,1])
        ax.set_yticks([0,1])                    
        ax.legend(['NWP_XXX','Guess_XXX','NWP_YYY','Guess_YYY'])
    fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.5)
    fig.suptitle(txt)