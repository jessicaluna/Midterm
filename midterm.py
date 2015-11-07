# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 19:55:14 2015

@author: jessicaluna
"""

from pylab import*
import numpy as np
import matplotlib.pyplot as plt
#from astropy.time import Time
#import csv
#import glob
import time
#import scipy.stats

starttime=time.time()

JD_time,RV_obs, RV_obs_err=loadtxt('RVmeasurements.txt',unpack=True,usecols=[0,1,2])
##function that calculates the sin curve for given parameters
#msini=ampllitude, and tzero is the offset
def rv_curve(Amplitude,period,tzero,x):
    RV_fit= Amplitude*np.sin((2*pi*(x-tzero)/period))
    return RV_fit
        
def newstep(amp,per,offs,ns):
    ha=4.#10.0**(0)
    hp=.01#10.0**(-1)
    ho=1.#10.0**(-1)
    a_step=amp
    p_step=per
    off_step=offs
    for kk in range(ns):    
        a_step=a_step +np.random.uniform(-ha,ha)
        p_step=p_step+ np.random.uniform(-hp,hp)
        off_step=off_step + np.random.uniform(-ho,ho)
    return(a_step,p_step,off_step)
    
def chisqfit(Rvobs,Rvmod,Rvmod_trial):
    std_current=np.std(RVmod)
    std_trial=np.std(Rvmod_trial)
    
    chi=sum(((Rvobs-Rvmod)/std_current)**2)
    chi_trial=sum(((Rvobs-Rvmod_trial)/std_trial)**2)
    print(" CHI:", chi,chi_trial)
    ratio=exp(-0.5*(chi_trial-chi))
    alph=min(1,ratio)
    return (alph)



    
    
    
def MCMCburners(Ab,Periodb,offsetb,niters,nburners):
    naccept=0.
    naccept_p=0.
    naccept_o=0.
    A_samples = []
    P_samples = []
    off_samples= []
    
    n_actual=niters
    for jj in range(niters+nburners):   
        RVcurrentb=rv_curve(Ab,Periodb,offsetb,JD_time)

        A_trialb, P_trialb, Offset_trialb=newstep(Ab,Periodb,offsetb,1)
        RV_trialb=rv_curve(A_trialb,P_trialb,Offset_trialb,JD_time)
        alphaa=chisqfit(RV_obs,RVcurrentb,RV_trialb)

        ua=np.random.uniform(0.0,1.0)
        if ua <=alphaa:
            Ab=A_trialb
 
            if nburners<jj:
                naccept +=1
                A_samples.append(Ab)
#        ## changing period        
        RVcurrentb=rv_curve(Ab,Periodb,offsetb,JD_time)
        RV_trialb=rv_curve(Ab,P_trialb,Offset_trialb,JD_time)
        alphab=chisqfit(RV_obs,RVcurrentb,RV_trialb)           

        ub=np.random.uniform()
        if ub <=alphab:
            Periodb=P_trialb
            if nburners<jj:
                naccept_p +=1
                P_samples.append(Periodb)
         ## changing zeropoint        
        RVcurrentb=rv_curve(Ab,Periodb,offsetb,JD_time)
        RV_trialb=rv_curve(Ab,Periodb,Offset_trialb,JD_time)
        
        alphac=chisqfit(RV_obs,RVcurrentb,RV_trialb)       

        uc=np.random.uniform()
        if uc <=alphac:
            offsetb=Offset_trialb 
            if nburners<jj:
                naccept_o +=1
                off_samples.append(offsetb)
                
    print('naccept',naccept,naccept_p,naccept_o, len(A_samples))
    print ('efficiency',float(naccept_o)/float(n_actual))
    print('last values',Ab,Periodb,offsetb)
    return(Ab,Periodb,offsetb,A_samples,P_samples,off_samples)


  
nburn=100
niter=10000

#initial guess
Astart=80.0
Periodstart=3.3
offsetstart=2452850.0

A,Period,offset,Asamples,Psamples,offsetsamples=MCMCburners(Astart,Periodstart,offsetstart,niter,nburn)
print (A,Period,offset)

savetxt('parameters.txt', zip(Asamples,Psamples,offsetsamples), fmt='%.18e', delimiter=' ', newline='\n', header='A,P,Tzero')



RVmodel=rv_curve(A,Period,offset,JD_time)

figure()
plt.errorbar(JD_time, RV_obs, yerr=RV_obs_err)
plot(JD_time,RVmodel,c='r',marker='o',ls='none')
xlabel('Time in JD')
ylabel('RV m/s')
title('RV vs Time')
show()

print(time.time()-starttime)
