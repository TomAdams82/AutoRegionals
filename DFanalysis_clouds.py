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

def simpleScoreDF(returnWorkings,TP_step,CC_step,HCC_step,guessDF_TP,guessDF_CC,guessDF_HCC,DF_priority,SDF,isWholeArea,inArea,TP_bins,CC_bins,XorY): # TP_step is TP at t
                    
    # Take a combination of guessed XXX or YYY and SSS and give it a score based on NWP -truth for that time
    
#    if ((XXXn == 0) & (YYYn == 5) & (TTTn == 1)):
#        pdb.set_trace()

    # Seperate NWP into subregions for XXX and YYY
    if isWholeArea == False: # i.e. mixed contribution - only calc for part of area
        subregion_YYY = SDF
        if XorY == 'Y':
            nwpDF_TP , b2 = np.histogram(TP_step[subregion_YYY], bins = TP_bins, density = False)/(np.sum(subregion_YYY).astype('float')) # do not normalise, just divide by number of points
            nwpDF_CC , b2 = np.histogram(CC_step[subregion_YYY], bins = CC_bins, density = False)/(np.sum(subregion_YYY).astype('float')) # do not normalise, just divide by number of points
            nwpDF_HCC , b2 = np.histogram(HCC_step[subregion_YYY], bins = CC_bins, density = False)/(np.sum(subregion_YYY).astype('float')) # do not normalise, just divide by number of points
        else: # i.e. is XXX
            subregion_XXX = inArea * np.logical_not(subregion_YYY) # Must be in main area, then inverse of subarea
            # Check still adds up to whole area
#            if np.sum(subregion_XXX+subregion_YYY) != np.sum(inArea):
#                print('Error - Subregions dont add up - check code')
#                sys.exit()
            nwpDF_TP , b2 = np.histogram(TP_step[subregion_XXX], bins = TP_bins, density = False)/(np.sum(subregion_XXX).astype('float')) # do not normalise, just divide by number of points
            nwpDF_CC , b2 = np.histogram(CC_step[subregion_XXX], bins = CC_bins, density = False)/(np.sum(subregion_XXX).astype('float')) # do not normalise, just divide by number of points
            nwpDF_HCC , b2 = np.histogram(HCC_step[subregion_XXX], bins = CC_bins, density = False)/(np.sum(subregion_XXX).astype('float')) # do not normalise, just divide by number of points
            if np.sum(subregion_XXX)==0:
                pdb.set_trace()
    else: # Whole area
        nwpDF_TP , b2 = np.histogram(TP_step[inArea], bins = TP_bins, density = False)/(np.sum(inArea).astype('float')) # do not normalise, just divide by number of points
        nwpDF_CC , b2 = np.histogram(CC_step[inArea], bins = CC_bins, density = False)/(np.sum(inArea).astype('float')) # do not normalise, just divide by number of points
        nwpDF_HCC , b2 = np.histogram(HCC_step[inArea], bins = CC_bins, density = False)/(np.sum(inArea).astype('float')) # do not normalise, just divide by number of points
    score_TP = KL_divergence(nwpDF_TP, guessDF_TP, base=None)       
    score_CC = KL_divergence(nwpDF_CC, guessDF_CC, base=None)       
    score_HCC = KL_divergence(nwpDF_HCC, guessDF_HCC, base=None)
    if DF_priority[2] == 0: # ignore high cloud, avoid 0 * inf = nan
        score = np.sum(DF_priority[:2] * np.array([score_TP,score_CC])) # scale scores by prioirty
    else:
        score = np.sum(DF_priority * np.array([score_TP,score_CC,score_HCC])) # scale scores by prioirty
    if score == np.nan:
        print('0 times infinity?')
        pdb.set_trace()
    if returnWorkings:
        return(score,score_TP,score_CC,score_HCC,nwpDF_TP,guessDF_TP,nwpDF_CC,guessDF_CC,nwpDF_HCC,guessDF_HCC)
    else:
        return(score)

def simpleScoreDF_winds(returnWorkings,U_step,V_step,guessDF_U,guessDF_V,DF_priority,SDF,isWholeArea,inArea,W_bins,XorY): # TP_step is TP at t
                    
    # Take a combination of guessed XXX or YYY and SSS and give it a score based on NWP -truth for that time - for winds
    
#    if ((XXXn == 0) & (YYYn == 5) & (TTTn == 1)):
#        pdb.set_trace()

    # Seperate NWP into subregions for XXX and YYY
    if isWholeArea == False: # i.e. mixed contribution - only calc for part of area
        subregion_YYY = SDF
        if XorY == 'Y':
            nwpDF_U , b2 = np.histogram(U_step[subregion_YYY], bins = W_bins, density = False)/(np.sum(subregion_YYY).astype('float')) # do not normalise, just divide by number of points
            nwpDF_V , b2 = np.histogram(V_step[subregion_YYY], bins = W_bins, density = False)/(np.sum(subregion_YYY).astype('float')) # do not normalise, just divide by number of points
        else: # i.e. is XXX
            subregion_XXX = inArea * np.logical_not(subregion_YYY) # Must be in main area, then inverse of subarea
            # Check still adds up to whole area
#            if np.sum(subregion_XXX+subregion_YYY) != np.sum(inArea):
#                print('Error - Subregions dont add up - check code')
#                sys.exit()
            nwpDF_U , b2 = np.histogram(U_step[subregion_XXX], bins = W_bins, density = False)/(np.sum(subregion_XXX).astype('float')) # do not normalise, just divide by number of points
            nwpDF_V , b2 = np.histogram(V_step[subregion_XXX], bins = W_bins, density = False)/(np.sum(subregion_XXX).astype('float')) # do not normalise, just divide by number of points
            if np.sum(subregion_XXX)==0:
                pdb.set_trace()
    else: # Whole area
        nwpDF_U , b2 = np.histogram(U_step[inArea], bins = W_bins, density = False)/(np.sum(inArea).astype('float')) # do not normalise, just divide by number of points
        nwpDF_V , b2 = np.histogram(V_step[inArea], bins = W_bins, density = False)/(np.sum(inArea).astype('float')) # do not normalise, just divide by number of points
    score_U = KL_divergence(nwpDF_U, guessDF_U, base=None)       
    score_V = KL_divergence(nwpDF_V, guessDF_V, base=None)       
    score = np.sum(DF_priority * np.array([score_U,score_V])) # scale scores by prioirty
    if returnWorkings:
        return(score,score_U,score_V,nwpDF_U,guessDF_U,nwpDF_V,guessDF_V)
    else:
        return(score)

    
def scoreDF(returnWorkings,t,TP,DF,TDF,TDF_text,SDF,SDF_text,XXXn,YYYn,TTTn,SSSn,inArea,TP_bins):  # OBSOLETE
    
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
        
def getText(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn, ruleCombos,ruleReplaces):
    
    initTxt = TDF_text[TTTn] # Initial text form before putting in weather
        
    if SDF_text[SSSn] == 'All': # i.e. a whole region change, no subregions
        txt = initTxt.replace('XXX',DF_text[XXXn]).replace('YYY',DF_text[YYYn])
    else:
        spatial_text = DF_text[YYYn] + ' ' + SDF_text[SSSn]
        txt = initTxt.replace('XXX',DF_text[XXXn]).replace('YYY',spatial_text)        
        # Make any changes accoridng to rules
    for i in range(len(ruleCombos)):
        if XXXn == ruleCombos[i,0]:
            if YYYn == ruleCombos[i,1]:
                if ruleCombos[i,2] == 2: #2nd occurence
                    txt = txt.replace(ruleReplaces[i,0],'$$$',1).replace(ruleReplaces[i,0],ruleReplaces[i,1],2).replace('$$$',ruleReplaces[i,0])                
                else: # 1st occurence
                    txt = txt.replace(ruleReplaces[i,0],ruleReplaces[i,1])
    # check that if have 'more frequent / widespread in the afternoon' change 'in the' to 'from'
    # In that case, replace TDF part with from
    if np.logical_or((np.char.find(txt,'widespread')>0),(np.char.find(txt,'more frequent')>0)): # check if contains widespread
        loc_inthe = np.char.find(TDF_text[TTTn] ,'in the')
        if loc_inthe > 0: # check TDF contains 'in the', but not 'and' e.g. 'in the afternoon and evening
            if np.char.find(TDF_text[TTTn] ,'and')<0: # check TDF contains 'in the', but not 'and' e.g. 'in the afternoon and evening
                search_txt = TDF_text[TTTn][loc_inthe:]
                new_txt = search_txt.replace('in the','from')
                txt = txt.replace(search_txt,new_txt)                       
        
    # Make first letter uppercase
    txtOut = txt[0].upper()+txt[1:]
    return txtOut


def getText_winds(DF_text,TDF_text,SDF_text,XXXn,YYYn,TTTn,SSSn):
    
    initTxt = TDF_text[TTTn] # Initial text form before putting in weather
    YYYtext = DF_text[YYYn]
    # Rules for winds (not from spreadsheet yet - too much repitition)
    strengthXXX = np.floor((XXXn+7)/8).astype(int)
    strengthYYY = np.floor((YYYn+7)/8).astype(int)
    drnXXX = XXXn - (strengthXXX-1)*8
    drnYYY = YYYn - (strengthYYY-1)*8
    if np.char.find(initTxt,'becoming YYY')>0: # check if contains
        if strengthXXX == 0: # change becoming to developing if light winds
            initTxt = initTxt.replace('becoming YYY','YYY developing')
        else:
            if strengthXXX < strengthYYY: # change becoming to developing if light winds
                if drnXXX == drnYYY: # Same direction
                    if strengthYYY<6:# Not Severe Gale, which has two spaces
                        YYYtext = YYYtext[:YYYtext.find(' ')] # drop anything after the space
                    else:
                        YYYtext = YYYtext[:(YYYtext.find('gale ')+4)] # drop anything before the direction
                else:
                    initTxt = initTxt.replace('becoming YYY','rising to YYY')
            if strengthXXX == strengthYYY: # change becoming to developing if light winds
                if drnXXX == drnYYY: # Same direction
                    print('Invalid change!')
                    pdb.set_trace()
                else:
                    initTxt = initTxt.replace('becoming YYY','turning YYY')
                    if strengthYYY<6:# Not Severe Gale, which has two spaces
                        YYYtext = YYYtext[(YYYtext.find(' ')+1):] # drop anything before the space
                    else:
                        YYYtext = YYYtext[(YYYtext.find('gale ')+1):] # drop anything before the direction
                    YYYtext = YYYtext.replace('lies','ly') # change lies to ly
            if strengthXXX > strengthYYY: # change becoming to developing if light winds
                if drnXXX == drnYYY: # Same direction
                    if strengthYYY == 1:
                        initTxt = initTxt.replace('becoming YYY','easing to a light breeze')    
                    else:
                        if strengthYYY == 2:
                            initTxt = initTxt.replace('becoming YYY','easing')    
                        else:
                            YYYtext = YYYtext[:YYYtext.find(' ')] # drop anything after the space
                else:
                    initTxt = initTxt.replace('becoming YYY','easing to YYY')                    


    if SDF_text[SSSn] == 'All': # i.e. a whole region change, no subregions
        txt = initTxt.replace('XXX',DF_text[XXXn]).replace('YYY',YYYtext)
    else:
        spatial_text = YYYtext + ' ' + SDF_text[SSSn]
        txt = initTxt.replace('XXX',DF_text[XXXn]).replace('YYY',spatial_text)      
    
    # Make first letter uppercase
    txtOut = txt[0].upper()+txt[1:]
    return txtOut

def getIcon(XXXmatrix,weightings,DF_icons): #input XXXmatrix is the score for when X was whole area , dims XXX type vs. t   
    totalScore = np.dot(XXXmatrix,weightings)
    best = np.where(totalScore == np.min(totalScore))[0]
    return DF_icons[best][0]    