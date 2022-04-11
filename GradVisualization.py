# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 09:56:54 2022

@author: lridings
"""

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

class GradVis:
    
    def __init__(self):
        self.something=0
    
    def cost_plot(self,results,cost_lw=1,legend=True,ax=None,ax2=None):
        if not ax:
            fig,ax=plt.subplots()
            ax2=ax.twinx()
        costHist=results[1]
        thetaHist=results[2]
        ax.plot(costHist,c='black',lw=cost_lw,label='Cost')
        ax.set_ylabel('Cost/RSS')
        ax2.set_ylabel('Theta/Weights')
        for el in list(thetaHist.columns):
            ax2.plot(thetaHist[el],ls='dotted',alpha=.4,label=el)
        if legend:
            ax.legend(loc='upper left')
            ax2.legend(loc='upper right')
        return ax,ax2
    def cost_summary(self,result_list):
        ax,ax2=self.cost_plot(result_list[0],legend=False)
        [self.cost_plot(result,cost_lw=.2,legend=False,ax=ax,ax2=ax2) for result in result_list]
    
    def theta_hist(self,result_list):
        theta_final = {}
        for el in result_list[0][0]:
            theta_final[el]=[]
        for result in result_list:
            for el in result[0]:
                theta_final[el]+=[result[0][el]]
        ax=sns.histplot(pd.DataFrame(theta_final))
        ax.set_xlabel('Theta/Weight')
        return theta_final
    def target_vs_guess(self,target,guess):
        fig,ax=plt.subplots()
        ax.plot(guess[:,0],target,label='Target')
        ax.plot(guess[:,0],guess[:,1],label='Guess')
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Intensity')
        ax.legend()
    
    def guess_resid_plot(self,target,guess):
        fig,ax=plt.subplots()
        resids=target-guess[:,1]
        ax.scatter(guess[:,0],resids)
        ax.plot(guess[:,0],guess[:,0]*0,color='red')
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Distance from true value')
        return fig,ax
    def plot_finalCosts(self,result_list):
        ax=sns.histplot(data=[r[1][-1] for r in result_list])
        return ax