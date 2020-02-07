import streamlit as st
# To make things easier later, we're also importing numpy and pandas for
# working with sample data.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import time 
import altair as alt
import scipy.special
from sidebar import run_sidebar


#%%






#%%
###%%% bin first then display
st.title('traveling fitness waves')
st.text('Matthew Melissa, PhD Candidate at Harvard University')
st.markdown('modeling the interplay of **mutation**, **selection** and **random genetic drift** in shaping the course of evolution') 

N,Ub,Ud,sb,sd,draw_b,draw_d,num_gen = run_sidebar()

def runsim(N,Ub,Ud,draw_b,draw_d,num_gen,seed):
	latest_iteration = st.empty()
	#lineages = [(1,)]
	extant_lineages = np.array([0]).astype(int)
	sizes = np.zeros((int(2*N*(Ub+Ud)*num_gen),num_gen))
	sizes[0,0] = N
	fits =  np.array([1.0])
	children = [[]]
	parents = [()]
	#big_lineages = np.array([]).astype(int)

	#progress_interval = num_gen/50
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
	            children.append([])
	            children[j] += [len(fits)-1]
	            parents.append(())
	            parents[-1] +=(j,)
	            sizes[len(fits)-1,t+1] = 1.0
	            sizes[j,t+1] -= 1.0
	    for j in extant_lineages:
	    	new_lineages_d = np.random.binomial(sizes[j,t+1],Ud)
	    	for k in range(new_lineages_d):
	            #lineages.append(lineages[j] + (np.random.randint(10**10),))
	            fits = np.append( fits,np.exp(-draw_d())*fits[j] )
	            extant_lineages = np.append(extant_lineages,len(fits))
	            children.append([])
	            children[j] += [len(fits)-1]
	            parents.append(())
	            parents[-1] +=(j,)
	            sizes[len(fits)-1,t+1] = 1.0
	            sizes[j,t+1] -= 1.0
	            
	            
	    extant_lineages = extant_lineages[sizes[extant_lineages,t+1]>0]

	    #big_lineages = np.union1d(big_lineages,extant_lineages[sizes[extant_lineages,t+1]>=10])


	    #if np.mod(t/2.0,progress_interval)==0 and t/progress_interval/2 <1:
	    	#bar.progress(t/progress_interval/2)


	#if len(big_lineages)==1:
		#big_lineages = np.arange(1,len(lineages))


	return children,parents,sizes[:len(children),:],fits

def add_trads(children,sizes,lineage_no):
    if len(children[lineage_no])==0:
        return sizes[lineage_no,:]
    else:
        sumsizess = np.copy(sizes[lineage_no,:])
        for ki in range(len(children[lineage_no])):
            sumsizess += add_trads(children,sizes,children[lineage_no][ki])
        return sumsizess



def add_trads_parents(children,parents,sizes):
	sumsizes = np.copy(sizes)
	no_children_pool = [i for i in range(len(children)) if len(children[i])==0]
	while True:
		for i in range(len(no_children_pool)):
			rem_individual = no_children_pool.pop(0)
			rem_parent = parents[rem_individual][0] 
			sumsizes[rem_parent,:] += sumsizes[rem_individual,:]
			children[rem_parent].remove(rem_individual)
			if len(children[rem_parent]) == 0:
				no_children_pool.append(rem_parent)

		if no_children_pool==[0] or len(no_children_pool)==0:
			break
	return sumsizes[:,:]






def plot_trajectories(N,num_gen,children,parents,sizes):
	chart_data = pd.DataFrame([])
	ki = 1
	#arrlist = [list(add_trads(children,sizes,k)) for k in big_lineages[1:]]
	#wide_df = pd.DataFrame(arrlist,columns=range(num_gen))


	sumsizes = add_trads_parents(children,parents,sizes)
	wide_df = pd.DataFrame(sumsizes,columns=range(num_gen))

	wide_df['lineage'] = np.arange(len(fits)).astype(str)

	chart_data = wide_df.melt(id_vars='lineage',var_name='generation',value_name='number individuals')
	chart_data=chart_data[chart_data['number individuals']>0]
	#chart_data = chart_data[chart_data['lineage']!='0']


	gen_range = pd.DataFrame({'generation':np.linspace(0,num_gen,100).astype(int)})
	c = alt.Chart(chart_data,width = 400).mark_line().encode(
		x='generation:Q',

		y=alt.Y('number individuals', scale=alt.Scale(domain=[0.0,1.1*N], clamp=True), title='number of individuals'),
		color=alt.Color('lineage',legend=None)
	)

	line = alt.Chart(chart_data).mark_line().encode(
		x='generation',

		y=alt.Y('number individuals', scale=alt.Scale(domain=[0.0,1.1*N], clamp=True), title='number of individuals'),
		color=alt.Color('lineage',legend=None)
	).interactive()

	nearest = alt.selection(type='single', nearest=True, on='mouseover',fields=['generation'], empty='none',init={'generation':0})


	#slider = alt.binding_range(min=0, max=int(num_gen)-1, step=10, name='generation:')
	#selector = alt.selection_single(name="generation", fields=['generation'],
                                #bind=slider, init={'generation': 50})

	# rules = alt.Chart(gen_range).mark_rule(color='gray').encode(
 #    	x='generation',
	# 	).transform_filter(selector)

	selectors = alt.Chart(gen_range).mark_point().encode(
    	x='generation',
    	opacity=alt.value(0),
		).add_selection(nearest)

	rules = alt.Chart(gen_range).mark_rule(color='gray').encode(
    	x='generation',
		).transform_filter(nearest)

	layerChart = alt.layer(line,selectors,rules).properties(width=600,height=300)

	df = pd.DataFrame(sizes[:,:],columns=range(num_gen))
	#df['fitness'] = df.index
	df['lineage'] = np.arange(len(fits)).astype(str)
	df['fitness'] = fits
	df_melt = df.melt(id_vars=['fitness','lineage'],var_name='generation',value_name='number')
	df_melt = df_melt[df_melt.number>0].astype({'generation': 'int32'})


	if Ub == 0:
		step_size = sd
	elif Ud == 0:
		step_size = sb
	else:
		step_size = np.min(sb,sd)

	bchart = alt.Chart(df_melt,width = 600).mark_bar().encode(
    	x=alt.Y('fitness:Q', bin=alt.Bin(extent=[min(df['fitness']), max(df['fitness'])], step=step_size), title='fitness'),
    	y=alt.Y('number:Q', scale=alt.Scale(domain=(0, N)), title='number of individuals'),
    	color=alt.Color('lineage',legend=None)).add_selection(nearest).transform_filter(nearest).interactive()


	concat_chart = (layerChart & bchart).configure_axis(labelFontSize=20,titleFontSize=20)
	return concat_chart



if st.sidebar.button('Run simulation'):
	seed = np.random.rand()
	bar = st.progress(0)
	#params = pd.DataFrame({'N':[N],'U':[U],'s':[s],'num_gen':[num_gen],'seed':[seed]})
	#params.to_csv('runparams.csv')


	#params = pd.read_csv('runparams.csv')
	#N,U,s,num_gen,seed = params['N'][0],params['U'][0],params['s'][0],params['num_gen'][0],params['seed'][0]

	children,parents,sizes,fits = runsim(N,Ub,Ud,draw_b,draw_d,num_gen,seed)
	bar.progress(50)

	layerChart = plot_trajectories(N,num_gen,children,parents,sizes)
	bar.progress(75)
	
	imagine_mutations = st.empty()

	st.altair_chart(layerChart  )
	imagine_measuring = st.empty()


	bar.progress(100)
	time.sleep(0.1)

	imagine_mutations.markdown('Imagine sequencing the DNA of all the individuals in a population (say, of *E. coli* or influenza viruses). The colored lines below show how many individuals carry each mutation over time.')

	imagine_measuring.markdown('One could also measure the reproductive ability, or fitness, of the population over time. Mouse over the upper plot to see how the population moves through "fitness space" as a *traveling wave*.')

	st.markdown('During my PhD, I have used a combination of analytical theory and simulations to describe the rate at which a population adapts and the shape of its traveling fitness wave. I have analyzed how these and other quantities depend on the size of the population, its mutation rate, and the distribution of fitness effects of new mutations, with a focus on large microbial populations.')

	st.markdown('Source code and more details are available at <https://github.com/mjmel/evo-vis>.')











