import numpy as np
import streamlit as st


@st.cache(hash_funcs={np.ufunc:id},suppress_st_warning=True)
def runsim(N,Ub,Ud,draw_b,draw_d,num_gen,assay_interval,seed):
	latest_iteration = st.empty()
	extant_lineages = np.array([0]).astype(int)
	sizes = np.zeros((int(2*N*(Ub+Ud)*num_gen),num_gen))
	sizes[0,0] = N
	fits =  np.array([1.0])
	children = [[]]
	parents = [()]

	assay_timepoints = np.arange(num_gen,step=assay_interval).astype(int)

	t=0
	for t in range(num_gen-1):
	    mean_fit = np.average(fits[extant_lineages],weights=sizes[extant_lineages,t])
	    sizes[extant_lineages,t+1] = np.random.poisson(sizes[extant_lineages,t]*np.exp(np.array(fits[extant_lineages])-mean_fit)*np.exp(1-1/float(N)*np.sum(sizes[extant_lineages,t])))
	    
	    for j in extant_lineages:
	        new_lineages = np.random.binomial(sizes[j,t+1],Ub)
	        for k in range(new_lineages):
	            fits = np.append( fits,np.exp(draw_b())*fits[j] )
	            extant_lineages = np.append(extant_lineages,len(fits))
	            children.append([])
	            children[j] += [len(fits)-1]
	            parents.append(())
	            parents[-1] +=(j,)
	            sizes[len(fits)-1,t+1] = 1.0
	            sizes[j,t+1] -= 1.0
	    for j in extant_lineages:
	    	new_lineages_d = np.random.binomial(sizes[j,t+1],Ud)
	    	for k in range(new_lineages_d):
	            fits = np.append( fits,np.exp(-draw_d())*fits[j] )
	            extant_lineages = np.append(extant_lineages,len(fits))
	            children.append([])
	            children[j] += [len(fits)-1]
	            parents.append(())
	            parents[-1] +=(j,)
	            sizes[len(fits)-1,t+1] = 1.0
	            sizes[j,t+1] -= 1.0
	                        
	    extant_lineages = extant_lineages[sizes[extant_lineages,t+1]>0]

	return children,parents,sizes[:len(children),assay_timepoints],fits,assay_timepoints