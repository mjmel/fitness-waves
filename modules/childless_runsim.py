import numpy as np
import streamlit as st

def runsim(N,Ub,Ud,draw_b,draw_d,num_gen,assay_interval,seed):
	latest_iteration = st.empty()
	#lineages = [(1,)]
	extant_lineages = np.array([0]).astype(int)
	sizes = np.zeros((int(2*N*(Ub+Ud)*num_gen),num_gen))
	sizes[0,0] = N
	fits =  np.array([1.0])
	num_children = []
	parents = []

	t=0
	for t in range(num_gen-1):
	    mean_fit = np.average(fits[extant_lineages],weights=sizes[extant_lineages,t])
	    sizes[extant_lineages,t+1] = np.random.poisson(sizes[extant_lineages,t]*np.exp(np.array(fits[extant_lineages])-mean_fit)*np.exp(1-1/float(N)*np.sum(sizes[extant_lineages,t])))
	    
	    for j in extant_lineages:
	        new_lineages = np.random.binomial(sizes[j,t+1],Ub)
	        for k in range(new_lineages):
	            #lineages.append(lineages[j] + (np.random.randint(10**10),))
	            fits = np.append( fits,np.exp(draw_b())*fits[j] )
	            extant_lineages = np.append(extant_lineages,len(fits))
	            num_children.append(0)
	            num_children[j] += 1
	            parents += [j]
	            sizes[len(fits)-1,t+1] = 1.0
	            sizes[j,t+1] -= 1.0
	    for j in extant_lineages:
	    	new_lineages_d = np.random.binomial(sizes[j,t+1],Ud)
	    	for k in range(new_lineages_d):
	            #lineages.append(lineages[j] + (np.random.randint(10**10),))
	            fits = np.append( fits,np.exp(-draw_d())*fits[j] )
	            extant_lineages = np.append(extant_lineages,len(fits))
	            num_children.append(0)
	            num_children[j] += 1
	            parents += [j]
	            sizes[len(fits)-1,t+1] = 1.0
	            sizes[j,t+1] -= 1.0
	            
	            
	    extant_lineages = extant_lineages[sizes[extant_lineages,t+1]>0]

	assay_timepoints = np.arange(num_gen,step=assay_interval).astype(int)
	assay_sizes = sizes[:len(num_children),assay_timepoints]



### pruning step. should probably run this every now and then
	to_keep = np.where((np.sum(assay_sizes,axis=1)>0) + (np.array(num_children)>0))[0].astype(int)
	num_kept = len(to_keep)
	keep_dict = {to_keep[i]:i for i in range(num_kept)}

	new_parents = np.zeros(len(to_keep)).astype(int)
	new_children = [[] for i in range(len(to_keep))]

	for i in range(len(to_keep)):
		new_parents[i] = keep_dict[parents[to_keep[i]]]
		new_children[new_parents[i]].append(i)

	parents = np.copy(new_parents)
	children = np.copy(new_children)

	assay_sizes = assay_sizes[to_keep,:]
	fits = fits[to_keep]


	return children,parents,assay_sizes,fits,assay_timepoints


