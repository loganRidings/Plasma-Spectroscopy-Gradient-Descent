# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 15:47:41 2022

@author: lridings
"""
#%% Imports
import LineID
import SpectralDataProcessing as sdp
import numpy as np
import pandas as pd
#%%
class WeightFinding:
    """
    Uses a gradient descent algorithm to find the best linear combination of characteristic
    spectra from various elements to fit a target spectrum. Control the species
    under consideration, the learning rate, and the iterations to get the best fit to a target spectrum.
    """
    def __init__(self):
        """
        Useful instance attributes
        Use spectral database via LineID
        Get the target spectrum from DataProcessing module
            reads in saved spectra from datatable
        center, std characterize the attenuation we apply to the database spectrum
            based on the measured attenuation of the experimental setup
        """
        self.pixelwidth = 0.253697
        self.numEl = 0
        self.data = sdp.DataProcessing()
        self.data.readIn()
        self.data.dropData()
        self.refLines = LineID.lineID()
        self.target = None
        self.setTargetSpectrum(norm=False)
        self.refLines.splitData()
        self.base = None
        self.attenuation = False
        self.center=700
        self.std=100
        self.elements=None
        self.pixels=None
        self.m_f=2
#%%
    def setBaseSpectrum(self):
        base = self.refLines.getBaseSpectrum(self.pixels,attenuate=self.attenuation,center=self.center,std=self.std,molecular_factor=self.m_f)
        self.base = base
    def getTargetSpectrum(self,wavelength):
        """Get the measured spectrum where 'wavelength' is at max intensity"""
        t = self.data.getPeakIntensity(wavelength)
        return self.data.getSpectrum(t)
    
    def setTargetSpectrum(self,wavelength=794.928,norm=True):
        """Set the target spectrum to the peaked intensity at wavelength
            bool norm will get normalized spectrum if true"""
        if norm: self.target = self.normalizeTarget(wavelength)
        else: self.target = self.getTargetSpectrum(wavelength)
        
    def normalizeTarget(self,wavelength=794.928,spec=None):
        """Get spectrum and normalize to area=1
            *beneficial to comparison between compositions where total measured brightness is different"""
        if spec is None:
            spec = self.getTargetSpectrum(wavelength)
        area = spec.sum()*self.pixelwidth
        spec = spec/area
        return spec
#%%
    def getGrad(self,species,ws):
        """Get the gradient for species given the total weighted and the base spectrums"""
        targ=self.target
        area = targ.sum()
        grad = -2*self.pixelwidth*np.sum((targ-ws[:,1])*self.base[species][:,1])/area
        return grad
    
    def cost(self,ws):
        """Cost function for the given weights"""
        area = self.target.sum()
        cost = self.pixelwidth*np.sum(np.square(ws[:,1]-self.target))/area
        return cost
    
    def gradientDescent(self,elements=['O I', 'O II','O III','O2','O2+'],learning_rate=1e-5,iterations=100):
        """
        Iteratively find the minimum cost, i.e., closest fit to the target spectrum via stepping
        weights to lowest cost configuration.
        default element list is all species of oxygen
        ------------------------
        Returns final weights, list of costs over iterations, and list of weights over iterations.
        """
        self.refLines.ellist=elements
        self.numEl=len(self.refLines.ellist)
        self.pixels=np.asarray(self.data.df.columns[3:].values,dtype='float')
        if not self.base:
            self.setBaseSpectrum()
        costHist = np.zeros(iterations)
        thetaHist=np.zeros((iterations,len(elements)))
        theta = {}
        for el in self.refLines.ellist:
            theta[el]=1
        for it in range(iterations):
            wspec = self.refLines.makeWeightedSpectrum(self.base, theta)
            for el in theta.keys():
                grad = self.getGrad(el,wspec)
                theta[el] = theta[el]-learning_rate*grad
            thetaHist[it,:]=list(theta.values())
            cost=self.cost(wspec)
            costHist[it]=cost
        return theta, costHist, pd.DataFrame(thetaHist,columns=elements)
