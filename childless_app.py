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
from childless_runsim import runsim

#%%

#%%
###%%% bin first then display
st.title('traveling fitness waves')
st.text('Matthew Melissa, PhD Candidate at Harvard University')
st.markdown('modeling the interplay of **mutation**, **selection** and **random genetic drift** in shaping the evolutionary dynamics of large populations') 

N,Ub,Ud,sb,sd,draw_b,draw_d,num_gen,assay_interval = run_sidebar()


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
			rem_parent = parents[rem_individual] 
			sumsizes[rem_parent,:] += sumsizes[rem_individual,:]
			children[rem_parent].remove(rem_individual)
			if len(children[rem_parent]) == 0:
				no_children_pool.append(rem_parent)

		if no_children_pool==[0] or len(no_children_pool)==0:
			break
	return sumsizes[:,:]

def plot_trajectories(N,num_gen,children,parents,sizes,assay_timepoints):
	chart_data = pd.DataFrame([])
	ki = 1
	sumsizes = add_trads_parents(children,parents,sizes)

	wide_df = pd.DataFrame(sumsizes,columns=assay_timepoints)
	wide_df['lineage'] = np.arange(len(fits)).astype(str)

	chart_data = wide_df.melt(id_vars='lineage',var_name='generation',value_name='number individuals')
	chart_data = chart_data[chart_data['number individuals']>0]

	last_zero = np.argmax(sumsizes!=0,axis=1)-1
	first_zero = sumsizes.shape[1] - (np.argmax(sumsizes[:,::-1]!=0,axis=1))
	plot_last_zero = (last_zero>=0)
	plot_first_zero = (first_zero<len(assay_timepoints))
	last_zero_df = pd.DataFrame({'lineage':np.arange(len(fits))[plot_last_zero].astype(str),'generation':assay_timepoints[last_zero[plot_last_zero]],'number individuals':np.zeros(len(fits))[plot_last_zero]})
	first_zero_df = pd.DataFrame({'lineage':np.arange(len(fits))[plot_first_zero].astype(str),'generation':assay_timepoints[first_zero[plot_first_zero]],'number individuals':np.zeros(len(fits))[plot_first_zero]})
	chart_data = chart_data.append(last_zero_df)
	chart_data = chart_data.append(first_zero_df)

	gen_range = pd.DataFrame({'generation':assay_timepoints})
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

	selectors = alt.Chart(gen_range).mark_point().encode(
    	x='generation',
    	opacity=alt.value(0),
		).add_selection(nearest)

	rules = alt.Chart(gen_range).mark_rule(color='gray').encode(
    	x='generation',
		).transform_filter(nearest)

	layerChart = alt.layer(line,selectors,rules).properties(width=600,height=300)

	df = pd.DataFrame(sizes[:,:],columns=assay_timepoints)
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
    	y=alt.Y('sum(number):Q', scale=alt.Scale(domain=(0, N)), title='number of individuals'),
    	color=alt.Color('lineage',legend=None)).add_selection(nearest).transform_filter(nearest).interactive()


	concat_chart = (layerChart & bchart).configure_axis(labelFontSize=20,titleFontSize=20)
	return concat_chart



if st.sidebar.button('Run simulation'):
	seed = np.random.rand()
	bar = st.progress(0)

	children,parents,sizes,fits,assay_timepoints = runsim(N,Ub,Ud,draw_b,draw_d,num_gen,assay_interval,seed)
	bar.progress(50)

	layerChart = plot_trajectories(N,num_gen,children,parents,sizes,assay_timepoints)
	bar.progress(75)
	
	imagine_mutations = st.empty()

	st.altair_chart(layerChart  )
	imagine_measuring = st.empty()


	bar.progress(100)
	time.sleep(0.1)

	imagine_mutations.markdown('Imagine sequencing the DNA of all the individuals in a population (say, of *E. coli* or influenza viruses). The colored lines below show how many individuals carry each mutation over time.')

	imagine_measuring.markdown('One can also measure the reproductive ability, or fitness, of the population over time. Mouse over the upper plot to see how the population moves through "fitness space" as a *traveling wave*.')

	st.markdown('During my PhD, I have used a combination of analytical theory and simulations to describe the rate at which a population adapts and the shape of its traveling fitness wave. I have analyzed how these and other quantities depend on the size of the population, its mutation rate, and the distribution of fitness effects of new mutations, with a focus on large microbial populations.')

	st.markdown('Source code and more details are available at <https://github.com/mjmel/evo-vis>.')











