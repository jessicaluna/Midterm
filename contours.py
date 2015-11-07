# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 13:47:54 2015

@author: jessicaluna
"""

from pylab import*
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#from astropy.time import Time
#import csv
#import glob
import time
import scipy.stats

starttime=time.time()

JD_time,RV_obs, RV_obs_err=loadtxt('RVmeasurements.txt',unpack=True,usecols=[0,1,2])

Aparams, Pparams, Offparams=loadtxt('parameters21chi8.txt',unpack=True,usecols=[0,1,2])

Msin_samples = msini(Aparams)
Period_samples=Pparams
Tzero_samples=Offparams

def rv_curve(Amplitude,period,tzero,x):
    RV_fit= Amplitude*np.sin((2*pi*(x-tzero)/period))
    return RV_fit

def msini(RVstar):
    Tday=3.5247 #days
    day_sec=86400 #sec in 1day
    T=Tday*day_sec
    Msun= 1.989*10**(33)
    Mstar=1.13* Msun
    G=6.67*10**(-8)
    
    R=((G*Mstar*T**2)/(4*pi**2))**(1.0/3.0)
    
    Vp=sqrt(G*Mstar/R)
    
    Mp= Mstar*RVstar/Vp
    return Mp

def stats(Sample):
    median=np.median(Sample)
    std=np.std(Sample)
    rangee=(max(Sample)-std)-(min(Sample)+std)
    return(median,std,rangee)    
    
#function compute_sigma_level taken from github
#http://jakevdp.github.io/blog/2014/06/14/frequentism-and-bayesianism-4-bayesian-in-python/    
def compute_sigma_level(trace1, trace2, nbins=20):
    """From a set of traces, bin by number of standard deviations"""
    L, xbins, ybins = np.histogram2d(trace1, trace2, nbins)
    L[L == 0] = 1E-16

    shape = L.shape
    L = L.ravel()

    # obtain the indices to sort and unsort the flattened array
    i_sort = np.argsort(L)[::-1]
    i_unsort = np.argsort(i_sort)

    L_cumsum = L[i_sort].cumsum()
    L_cumsum /= L_cumsum[-1]
    
    xbins = 0.5 * (xbins[1:] + xbins[:-1])
    ybins = 0.5 * (ybins[1:] + ybins[:-1])

    return xbins, ybins, L_cumsum[i_unsort].reshape(shape)
    

Msin_median,Msin_std,Msin_range =stats(Msin_samples)
print('Msini med,std,range',Msin_median/(1.898*10**(30)),Msin_std,Msin_range/(1.898*10**(30)))
Period_median,Period_std,Period_range=stats(Period_samples)
print('Period med,std,range',Period_median,Period_std,Period_range)
Tzero_median,Tzero_std,Period_range=stats(Tzero_samples)
print('Tzero med,std,range',Tzero_median,Tzero_std,Period_range)

#Msini_pdf=mlab.normpdf(Msin_samples,Msin_median,Msin_std)
#Period_pdf=mlab.normpdf(Period_samples,Period_median,Period_std)
#Tzero_pdf=mlab.normpdf(Tzero_samples,Tzero_median,Tzero_std)

Msin,Period,sig_MP=compute_sigma_level(Msin_samples, Period_samples)
Msin,Tzero,sig_MT=compute_sigma_level(Msin_samples, Tzero_samples)
Tzero,Period,sig_TP=compute_sigma_level(Tzero_samples, Period_samples)

Abest=77.47156475290217
Pbest= 3.5244657203569587
Obest=2452888.2328501563

Msintru=msini(Abest)
print('Msini in solar mass',Msintru/(1.898*10**(30)))

t=linspace(JD_time[0],JD_time[JD_time.size-1], 50000)
RVfit=rv_curve(Abest,Pbest,Obest,t)
test=rv_curve(Abest,Pbest,Obest,JD_time)

figure(1)
plt.contour( Msin,Period,transpose( sig_MP), levels=[0.68, 0.95,0.997])
ylabel('Period in JD ')
xlabel('Msini in g')

figure(2)
plt.contour( Msin,Tzero,transpose( sig_MT), levels=[0.68, 0.95,0.997])
ylabel('Tzero in JD ')
xlabel('Msini in g')

figure(3)
plt.contour( Tzero,Period,transpose( sig_TP), levels=[0.68, 0.95,0.997])
xlabel('Tzero in JD ')
ylabel('Period in JD')

figure(4)
plt.errorbar(JD_time, RV_obs, yerr=RV_obs_err,ls='none')
plot(t,RVfit,c='k')#,marker='o',ls='none')
plt.axis([2452848,2452848+10,-100,100])
xlabel('Time in JD')
ylabel('RV m/s')
title('RV vs Time')

figure(5)
plt.errorbar(JD_time, RV_obs, yerr=RV_obs_err)
plot(JD_time,test,c='c',marker='o',ls='none')

show()

print(time.time()-starttime)
