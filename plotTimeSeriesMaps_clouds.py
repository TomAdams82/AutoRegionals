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

def plotTimeSeriesMaps(param,data_timestep,lats,lons,NWP,sf,maplims,inArea,TDF,SDF_text,SDF,SSSn,txt): # NWP is TP, CC, or HCC
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
        if param == 0: #TP
            m.pcolormesh(x, y, NWP[:,:,t], vmin=0, vmax=6, cmap='Blues')#
            m.colorbar()
        if param == 1: #CC
            m.pcolormesh(x, y, NWP[:,:,t], vmin=0, vmax=1, cmap='Greys')#
        if param == 2: # HCC
            m.pcolormesh(x, y, NWP[:,:,t], vmin=0, vmax=1, cmap='Greens')#            
        if param == 3: # U
            m.pcolormesh(x, y, NWP[:,:,t], vmin=-30, vmax=30, cmap='coolwarm')#            
        if param == 4: # V
            m.pcolormesh(x, y, NWP[:,:,t], vmin=-30, vmax=30, cmap='coolwarm')#            

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

def plotTimeSeriesDFs(param,bins,nwpDF_XXX,nwpDF_YYY,guessDF_XXX,guessDF_YYY,scoreXXX,scoreYYY,data_timestep,txt):  
    fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
    gs = mpl.gridspec.GridSpec(3, 3)
    nlabels = len(bins)-1  
    param_labels = np.zeros(nlabels).astype(str)
    #param_labels[0]=''
    #param_labels[-1]=''
    utc_labels = np.concatenate([(np.arange(15,24,3)),(np.arange(0,15,3))]).astype(str)
    for i in range(nlabels):
        param_labels[i] = str(bins[i])+'-'+str(bins[i+1])
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
        ax.set_xticks(x)
        ax.set_xlim([-1,nlabels+1])
        if param == 0:
            ax.set_xlabel('Precip amount / mm')
        if param == 1:     
            ax.set_xlabel('Cloud fraction')
        if param == 2:
            ax.set_xlabel('High cloud fraction')
        if param == 3:
            ax.set_xlabel('U component / kt')
        if param == 4:
            ax.set_xlabel('V component / kt')

        ax.set_xticklabels(param_labels, rotation = 'vertical',fontsize=6)
        ax.set_ylim([0,1])
        ax.set_yticks([0,1])  
        
        ax.legend(['NWP_XXX','Guess_XXX','NWP_YYY','Guess_YYY'])
    fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.5)
    fig.suptitle(txt)
    
def compareTimeSeriesDFs(param,bins,nwpDF_XXX,guessDF_XXX,scoreXXX,DF_txt,bestOption,txt):  
    fig = plt.figure(figsize=(20, 20), facecolor='w', edgecolor='w')
    nx = int(np.ceil(np.sqrt(len(DF_txt))))
    ny = int(np.ceil(len(DF_txt)/(1.0*nx)))
    offset = (nx*ny)-len(DF_txt)
    gs = mpl.gridspec.GridSpec(ny, nx)
    nlabels = len(bins)-1  
    param_labels = np.zeros(nlabels).astype(str)
    for i in range(nlabels):
        param_labels[i] = str(bins[i])+'-'+str(bins[i+1])
    x = np.arange(nlabels)
    for o in range(len(DF_txt)):
        #score[i] = KL_divergence(nwpDF[:,t], guessDF[:,t], base=None) 
        ax=plt.subplot(gs[o+offset])    
        ax.plot(x,nwpDF_XXX,'b-')
        ax.plot(x,guessDF_XXX[o,:],'b:')
        if o == bestOption:
            ax.set_title(DF_txt[o]+' '+str(scoreXXX[o]),fontsize=10,fontweight='bold')
        else:
            ax.set_title(DF_txt[o]+' '+str(scoreXXX[o]),fontsize=10)
        #ax.set_yticks(x)
        ax.set_xticks(x)
        ax.set_xlim([-1,nlabels+1])
        if param == 0:
            ax.set_xlabel('Precip amount / mm')
        if param == 1:     
            ax.set_xlabel('Cloud fraction')
        if param == 2:
            ax.set_xlabel('High cloud fraction')
#        if param == 3:
#            ax.set_xlabel('U component / kt')
#        if param == 4:
#            ax.set_xlabel('V component / kt')

        ax.set_xticklabels(param_labels, rotation = 'vertical',fontsize=6)
        ax.set_ylim([0,1])
        ax.set_yticks([0,1])  
    fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.6)
    fig.suptitle(txt)