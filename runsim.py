import numpy as np
import streamlit as st

def runsim(N,U,s,num_gen):

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

	    big_lineages = np.union1d(big_lineages,extant_lineages[sizes[extant_lineages,t+1]>=10])

	if len(big_lineages)==1:
		big_lineages = np.arange(1,len(lineages))


	return big_lineages,children,sizes,fits

