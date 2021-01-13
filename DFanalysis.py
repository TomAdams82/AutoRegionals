#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 18:46:23 2021

@author: tadams@met.co.nz
"""

import numpy as np
from scipy.stats import entropy as KL_divergence
import pdb
import sys

def simpleScoreDF(returnWorkings,TP_step,guessDF,SDF,isWholeArea,inArea,TP_bins,XorY): # TP_step is TP at t
    
    # Take a combination of guessed XXX or YYY and SSS and give it a score based on NWP -truth for that time
    
#    if ((XXXn == 0) & (YYYn == 5) & (TTTn == 1)):
#        pdb.set_trace()

    # Seperate NWP into subregions for XXX and YYY
    if isWholeArea == False: # i.e. mixed contribution - only calc for part of area
        subregion_YYY = SDF
        if XorY == 'Y':
            nwpDF , b2 = np.histogram(TP_step[subregion_YYY], bins = TP_bins, density = False)/(np.sum(subregion_YYY).astype('float')) # do not normalise, just divide by number of points
        else: # i.e. is XXX
            subregion_XXX = inArea * np.logical_not(subregion_YYY) # Must be in main area, then inverse of subarea
            # Check still adds up to whole area
#            if np.sum(subregion_XXX+subregion_YYY) != np.sum(inArea):
#                print('Error - Subregions dont add up - check code')
#                sys.exit()
            nwpDF , b2 = np.histogram(TP_step[subregion_XXX], bins = TP_bins, density = False)/(np.sum(subregion_XXX).astype('float')) # do not normalise, just divide by number of points
            if np.sum(subregion_XXX)==0:
                pdb.set_trace()
    else: # Whole area
        nwpDF , b2 = np.histogram(TP_step[inArea], bins = TP_bins, density = False)/(np.sum(inArea).astype('float')) # do not normalise, just divide by number of points

    score = KL_divergence(nwpDF, guessDF, base=None)       
    if returnWorkings:
        return(score,nwpDF,guessDF)
    else:
        return(score)
    
def scoreDF(returnWorkings,t,TP,DF,TDF,TDF_text,SDF,SDF_text,XXXn,YYYn,TTTn,SSSn,inArea,TP_bins):
    
    # Take a combination of guessed XXX, YYY, TTT and SSS and give it a score based on NWP -truth for that time
    
#    if ((XXXn == 0) & (YYYn == 5) & (TTTn == 1)):
#        pdb.set_trace()

    scoreXXX, scoreYYY = [0,0]
    nwpDF_XXX,nwpDF_YYY = [0,0]
    guessDF_XXX,guessDF_YYY = [0,0]
    TP_step = TP[:,:,t] 
    # Seperate NWP into relevance values 0 - 1 for subregions for XXX and YYY
    if TDF[TTTn,t]==True: # i.e. no contribution from YYY
        subregion_XXX = inArea
        nwpDF_XXX , b2 = np.histogram(TP_step[subregion_XXX], bins = TP_bins, density = False)/(np.sum(subregion_XXX).astype('float'))
        guessDF_XXX = DF[XXXn,:,0]
        scoreXXX = KL_divergence(nwpDF_XXX, guessDF_XXX, base=None)
        score = scoreXXX
    else:
        subregion_YYY = SDF[:,:,SSSn]
        nwpDF_YYY , b2 = np.histogram(TP_step[subregion_YYY], bins = TP_bins, density = False)/(np.sum(subregion_YYY).astype('float')) # do not normalise, just divide by number of points
        guessDF_YYY = DF[YYYn,:,0]
        scoreYYY = KL_divergence(nwpDF_YYY, guessDF_YYY, base=None)
        
        if SDF_text[SSSn] != 'All': # i.e. add contribution from XXX
            subregion_XXX = inArea * np.logical_not(subregion_YYY) # Must be in main area, then inverse of subarea
            # Check still adds up to whole area
            if np.sum(subregion_XXX+subregion_YYY) != np.sum(inArea):
                print('Error - Subregions dont add up - check code')
                sys.exit()
            nwpDF_XXX , b2 = np.histogram(TP_step[subregion_XXX], bins = TP_bins, density = False)/(np.sum(subregion_XXX).astype('float')) # do not normalise, just divide by number of points
            # Create guess DF_XXX
            guessDF_XXX = DF[XXXn,:,0]
            scoreXXX = KL_divergence(nwpDF_XXX, guessDF_XXX, base=None) 
            score = np.mean([scoreXXX,scoreYYY])
        else:
            if TDF_text[TTTn] == 'XXX but YYY': # prevent Rain but showers with no spatial qualifier type scenarios
                score = np.inf
            else:
                score = scoreYYY 
        
    if returnWorkings:
        return(score,scoreXXX,scoreYYY,nwpDF_XXX,guessDF_XXX,nwpDF_YYY,guessDF_YYY)
    else:
        return(score)
        
def getText(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn):
    
    initTxt = TDF_text[TTTn] # Initial text form before putting in weather
    if ((XXXn < 5) & (YYYn >4)): # change becoming to developing if wx getting significantly worse
        if np.char.find(initTxt,'becoming YYY')>0: # check if countains it
            initTxt = initTxt.replace('becoming YYY','YYY developing')
        
    if SDF_text[SSSn] == 'All': # i.e. a whole region change, no subregions
        txt = initTxt.replace('XXX',DF_text[XXXn]).replace('YYY',DF_text[YYYn])
    else:
        spatial_text = DF_text[YYYn] + ' ' + SDF_text[SSSn]
        txt = initTxt.replace('XXX',DF_text[XXXn]).replace('YYY',spatial_text)        
    # Special rules for repeated words

    if ((DF_text[XXXn] == 'drizzle') & (DF_text[YYYn] == 'patchy drizzle')):
        txt = txt.replace('patchy drizzle','patchy')
    if ((DF_text[XXXn] == 'isolated showers') & (DF_text[YYYn] == 'showers')): # Replace second occurence of showers
        txt = txt.replace('showers','$$$',1).replace('showers','widespread',2).replace('$$$','showers')        
    if ((DF_text[XXXn] == 'a few showers') & (DF_text[YYYn] == 'showers')): # Replace second occurence of showers
        txt = txt.replace('showers','$$$',1).replace('showers','more frequent',2).replace('$$$','showers')
    if ((DF_text[XXXn] == 'isolated showers') & (DF_text[YYYn] == 'showers, some heavy')):
        txt = txt.replace('showers, some heavy','widespread and some heavy')
    if ((DF_text[XXXn] == 'a few showers') & (DF_text[YYYn] == 'showers, some heavy')):
        txt = txt.replace('showers, some heavy','more frequent and some heavy')
    if ((DF_text[XXXn] == 'showers') & (DF_text[YYYn] == 'showers, some heavy')):
        txt = txt.replace('showers, some heavy','some heavy')
    if ((DF_text[XXXn] == 'rain') & (DF_text[YYYn] == 'rain with heavy falls')):
        txt = txt.replace('rain with heavy falls','with heavy falls')
    return txt

def DF_banned_list(Pn,XXXn,YYYn,banned_list_precip):
    # Return whether combination of parameters in permissible
    allowed = True
    if XXXn == YYYn:
        allowed = False
    else:
        if Pn == 0: # Precip
            for b in banned_list_precip[XXXn]:
                if YYYn == b:
                    allowed = False
    return allowed
    