#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 22:25:26 2020

@author: matthewmelissa
"""
import numpy as np
import matplotlib.pyplot as plt
#%%

num_gen = 2000

N = 10**3
U = 0.00005
s = 0.02
n_arr = np.zeros((N,num_gen))

#%%

lineages = [(1,)]
sizes = [[N]] 
fits =  [1.0]
children = [()]
extant_lineages = [1]

mean_fit = 1.0
t=0
for t in range(num_gen-1):
    for j in extant_lineages:
        new_lineages = np.random.binomial(sizes[j,-1],U)
        for k in range(new_lineages):
            lineages.append(lineages[j] + (np.random.randint(10**10),))
            children.append(())
            children[j] += (len(lineages)-1,)
            fits.append( np.exp(np.random.exponential(s))*fits[j] )
            sizes.append([1.0])
            extant_lineages.append(len(lineages)-1)

            if sizes[j][-1]-1 > 0:
                sizes[j].append(sizes[j][-1]-1)
            else:
                extant_lineages.remove(j)
            
            

    mean_fit = np.average(fits,weights=sizes[:len(fits),t])
    sizes[:len(fits),t+1] = np.random.poisson(sizes[:len(fits),t]*np.exp(np.array(fits)-mean_fit)*np.exp(1-1/float(N)*np.sum(sizes[:,t])))



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



