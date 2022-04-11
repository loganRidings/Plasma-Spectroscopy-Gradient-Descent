# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 07:59:31 2022
Use this script to translate the raw text which can be output from the SpecLine
software into an excel sheet.
The input file should be a .txt file. The output data
can be written any specified directory (local or network)
@author: Logan Ridings
"""

import pandas as pd
import numpy as np
import os

class lineID:
    """
    Create an object which has methods to collect spectral lines from a database,
    format and write them to csv, and perform operations to prepare them for 
    fitting a measured spectrum.
    """
    pixelCount=3648
    def __init__(self):
        """Constructor method. Default is to pull lines from the .txt file
        and use the specified list of element species."""
        self.lines = self.dothething()
        # self.ellist = ['O I', 'O II','O III','O2', 'O2+']
        self.ellist = ['O I', 'O II','O III','O2','O2+']
        self.allLines=None
    def dothething(self, writefile = False):
        """Read the text file in and make the dataframe to store and access available
        lines from."""
        # pathin = input("Filepath including extension: ")
        pathin = r"\\opdata2\Logan\Plasma Diagnostics\References\OFeNHlines.txt"
        pathin = pathin.replace("\"","")
        lines = []
        with open(pathin) as file:
            lines = file.readlines()
        cleanlines = []
        curWave = 0
        curInt = 0
        curLine = []
        for line in lines[8:]: # iterate through all lines from file
            if len(line) <=1: continue
            splitL = line.split("\t")
            for part in splitL:
                if len(part) == 0: 
                    continue
                if part.startswith("Wavelength:"): # update current measured wavelength
                    curWave = float(splitL[splitL.index(part)+1][:-4])
                    break
                if part.startswith("Intensity:"): # updata current measured intensity
                    curInt = float(splitL[splitL.index(part)+2])
                    break
                if part.startswith("-"):
                    break
                if part.startswith("Lines"):
                    break
                if part.startswith("Element"):
                    break
                if part.startswith('lower'):
                    break
                if 'lower' in splitL[4]: break
                # If the loop has gotten here, the current line is a spectral line
                # Add the important information to the list of measured lines
                element = splitL[2]
                nomLine = splitL[1]
                relInt = splitL[3]
                lowerE = splitL[4]
                upperE = splitL[6]
                lowerT=splitL[7]
                upperT=splitL[9]
                quantumN=splitL[10]
                comment = ''
                if len(splitL) ==14 and splitL: 
                    comment = splitL[-1]
                curLine = [curWave,curInt,element,nomLine,relInt,lowerE,upperE,lowerT,upperT,quantumN,comment]
                cleanlines.append(curLine)
                break
        df = pd.DataFrame(cleanlines,columns=['Measured Line [nm]','Intensity [counts]','Element',
                                              'Line [nm]','Rel Intensity','Lower Energy [eV]','Upper Energy [eV]',
                                              'lower state','upper state','quantum num','Comment'])
        df['Line [nm]']=df['Line [nm]'].astype(float)
        df['Rel Intensity']=df['Rel Intensity'].astype(float)
        df['Element']=df['Element'].str.strip()
        if writefile:
            writepath = input('Directory to save Data: ')
            writepath = r'\\opdata2\Logan\Plasma Diagnostics\References'
            filename = input("Filename for data (no extension): ")
            filename = 'ref'
            os.chdir(writepath)
            df.to_excel(filename+'.xlsx', index=False)
        return df
    def set_elements(self,elements: list):
        self.ellist=elements
        
    def getElement(self,element):
        """Return all species lines associated with the given element"""
        if element=='o': return self.lines[self.lines['Element'].isin(['O I', 'O II','O III', 'O2', 'O2+'])]
        if element=='fe': return self.lines[self.lines['Element'].isin(['Fe II'])]
        if element=='n': return self.lines[self.lines['Element'].isin(['N III'])]
    
    def splitData(self):
        """Get the subset of spectral lines from the database of interest"""
        o = self.getElement('o')
        fe = self.getElement('fe')
        n = self.getElement('n')
        elements = pd.concat([o,fe,n])
        elements.sort_values('Line [nm]',inplace=True)
        self.allLines = elements
    
    def getPixelLines(self,pixel,molecular_factor):
        """
        Turn each set of species lines into a spectrum on the domain of the pixel
        wavelengths which are measured on the spectrometer. The intensity of the species
        at pixel 'p' is given by the sum of the intensities of all species lines,
        the weight of each line determined by a how far the line is from the wavelength
        of the pixel (using a gaussian estimator). This weight is multiplied by the
        relative intensity of the line. The fwhm is representative of the resolution
        of the spectrometer and allows for a line of pixel wavelength+-1nm to be counted towards the
        intensity at pixel p
        ---------------------------------------
        Returns a dict of each species with an associated intensities
            for an input of 'pixel'
        """
        intensities = {}
        lines = self.allLines
        fwhm = 1
        sigma=fwhm
        for el in self.ellist:
            elInt = lines[lines['Element']==el]['Rel Intensity'].to_numpy()
            elLines = lines[lines['Element']==el]['Line [nm]'].to_numpy()
            if el in ['O2','O2+']:
                elInt = np.sum(elInt*np.exp(-np.square(elLines-pixel)/(2*(molecular_factor*sigma)**2))) #wide molecular bands
            else: elInt = np.sum(elInt*np.exp(-np.square(elLines-pixel)/(2*sigma**2))) #narrow atomic bands
            intensities[el]=elInt
        return intensities

    def getBaseSpectrum(self,pixels,attenuate=True,center=600,std=150,molecular_factor=1):
        """
        Use 'getPixelLines' on each of the pixels/wavelengths passed into the func
        and return a dict with species as keys and arrays with pixels,intensities
        as values
        If attenuate is true, pass to attenuation func.
        """
        spec = {}
        
        for i in range(len(pixels)):
            pix = self.getPixelLines(pixels[i],molecular_factor)
            for el in pix.keys():
                if not el in spec: spec[el]=[]
                spec[el].append(pix[el])
        for el in spec.keys():
            spec[el]=np.array([pixels,spec[el]]).T
        if attenuate:
            for el in spec.keys():
                spec[el]=self.attenuateSpectrum(spec[el],center,std)
        return spec
    
    def getArea(self,specArray):
        """Based on a given spectrum, return the area under the curve. Uses pixel width of spectrometer"""
        area = 0.253697*specArray[:,1].sum()
        return area
    
    def makeWeightedSpectrum(self,specdict,weights):
        """Multiply each species spectrum by the species weight. Return summed weighted spectrum."""
        wspec = np.zeros((self.pixelCount,2))
        for el in weights.keys():
            wspec[:,0]=specdict[el][:,0]
            wspec[:,1]+=specdict[el][:,1]*weights[el]
        return wspec
    
    def attenuateSpectrum(self,spec,center=700,std=100):
        """
        Shape spectrum based on attenuation present in optical path from plasma 
        through to spectrometer pixels.
        """
        var = std**2
        x = 1119.535*np.linspace(0,1,self.pixelCount)+194.68
        attenuation = np.exp(-np.square(x-center)/(2*var))
        spec[:,1]=spec[:,1]*attenuation
        return spec

def main():
    l = lineID()
    return l

if __name__=='__main__':
    l =main()