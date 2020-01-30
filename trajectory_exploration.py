#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 22:25:26 2020

@author: matthewmelissa
"""
import numpy as np
import matplotlib.pyplot as plt
%matplotlib osx
#%%

num_gen = 1000

N = 10**3
U = 0.001
s = 0.01



lineages = [(1,)]
extant_lineages = np.array([0]).astype(int)
sizes = np.zeros((int(2*N*U*num_gen),num_gen))
sizes[0,0] = N
fits =  np.array([1.0])
children = [()]
big_lineages = np.array([]).astype(int)


t=0
for t in range(num_gen-1):
    mean_fit = np.average(fits[extant_lineages],weights=sizes[extant_lineages,t])
    sizes[extant_lineages,t+1] = np.random.poisson(sizes[extant_lineages,t]*np.exp(np.array(fits[extant_lineages])-mean_fit)*np.exp(1-1/float(N)*np.sum(sizes[extant_lineages,t])))
    for j in extant_lineages:
        new_lineages = np.random.binomial(sizes[j,t+1],U)
        for k in range(new_lineages):
            lineages.append(lineages[j] + (np.random.randint(10**10),))
            extant_lineages = np.append(extant_lineages,len(lineages))
            children.append(())
            children[j] += (len(lineages)-1,)
            fits = np.append( fits,np.exp(np.random.exponential(s))*fits[j] )
            sizes[len(lineages)-1,t+1] = 1.0
            sizes[j,t+1] -= 1.0
            
    extant_lineages = extant_lineages[sizes[extant_lineages,t+1]>0]

    big_lineages = np.union1d(big_lineages,extant_lineages[sizes[extant_lineages,t+1]>25])


            
plt.figure()
for k in big_lineages[1:]:
    plt.plot(add_trads(k))
    
    
#%%
plt.figure()
for k in range(1,len(lineages)):
    plt.plot(add_trads(k))
    
#%%    

def add_trads(lineage_no):
    if len(children[lineage_no])==0:
        return sizes[lineage_no,:]
    else:
        sumsizess = np.copy(sizes[lineage_no,:])
        for ki in range(len(children[lineage_no])):
            sumsizess += add_trads(children[lineage_no][ki])
        return sumsizess

#%%    
    #%%
plt.figure()
upper = np.sum(sizes,axis=0)
lower = 0 
rgb = tuple(np.random.rand(3))

new_upper = upper - sizes[0,:]/2.0
new_lower = lower + sizes[0,:]/2.0
plt.fill_between(range(num_gen),new_lower,lower,color=rgb)
plt.fill_between(range(num_gen),new_upper,upper,color=rgb)

#%%

cumsizes = np.zeros((int(2*N*U*num_gen),num_gen))
for i in range(len(lineages)):
    for k in range(len(lineages)):
        
    cumsizes[i,:] +=



