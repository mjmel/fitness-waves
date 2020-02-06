import streamlit as st
# To make things easier later, we're also importing numpy and pandas for
# working with sample data.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import time 
import altair as alt
import scipy.special
#%%



@st.cache
def get_cache_seed():
	params = pd.DataFrame({'N':[100],'U':[0.01],'s':[0.01],'num_gen':[50],'seed':[1]})
	params.to_csv('runparams.csv')

	try:


		return cache_seed()

	except NameError:
		return np.random.rand()

st.text(get_cache_seed())

@st.cache()
def cache_seed():
	return 4



seed = 1


#%%
###%%% bin first then display
st.title('interactive evolution simulations')
st.text('Matthew Melissa, PhD Candidate at Harvard University')


st.markdown('evolution involves a complex interplay of **mutation**, **selection** and **genetic drift**') 

N = 10**st.sidebar.slider('choose a population size',min_value=2.0,max_value=4.0,step=0.01,value=3.2,format='10E%f')
U = 10**st.sidebar.slider('choose a mutation rate',min_value=-5.0,max_value=-2.0,step=0.01,value=-3.0,format='10E%f')


dfe_option = st.sidebar.selectbox('choose a distribution of fitness effects of mutations',pd.Series(['single effect','exponential','gamma']))

if dfe_option == 'gamma':

	gamma_shape_opt = st.sidebar.slider('choose a gamma shape parameter',min_value=1.0,max_value=32.0,step=1.0,value=1.0)



s = 10**st.sidebar.slider('choose an average effect size',min_value=-4.0,max_value=-1.0,step=0.01,value=-2.0,format='10E%f')

num_gen = st.sidebar.slider('how many generations?',min_value=20,max_value=10000,step=10,value=1000)


svals = np.linspace(0,4*s,100)

if dfe_option == 'exponential':
	dfe_source = pd.DataFrame({
	  'fitness effect': svals,
	  'probability': 1/s*np.exp(-svals/s)
	})
	def dfe():
		return np.random.exponential(s)


	c = alt.Chart(dfe_source,width=300).mark_area().encode(
    x='fitness effect',
    y='probability'

	)
	st.sidebar.altair_chart(c)

elif dfe_option == 'gamma':
	beta = gamma_shape/effect_size
	dfe_source = pd.DataFrame({
	  'fitness effect': svals,
	  'probability': beta**gamma_shape/scipy.special.gamma(gamma_shape)*svals**(gamma_shape-1)*np.exp(-beta*svals)
	})

	c = alt.Chart(dfe_source,width=300).mark_area().encode(
    x='fitness effect',
    y='probability'

	)
	st.sidebar.altair_chart(c)


elif dfe_option == 'single effect':
	dfe_source = pd.DataFrame({
	  'fitness effect': [s,s],
	  'probability': [0,1]
	})
	def dfe():
		return s
#st.sidebar.area_chart(dfe_source).encode(x='fitness effect',y='probability')

@st.cache(hash_funcs={np.ufunc:id},suppress_st_warning=True)
def runsim(N,U,dfe,num_gen,seed):
	latest_iteration = st.empty()
	lineages = [(1,)]
	extant_lineages = np.array([0]).astype(int)
	sizes = np.zeros((int(2*N*U*num_gen),num_gen))
	sizes[0,0] = N
	fits =  np.array([1.0])
	children = [()]
	big_lineages = np.array([]).astype(int)

	progress_interval = num_gen/50
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
	            fits = np.append( fits,np.exp(dfe())*fits[j] )
	            sizes[len(lineages)-1,t+1] = 1.0
	            sizes[j,t+1] -= 1.0
	            
	    extant_lineages = extant_lineages[sizes[extant_lineages,t+1]>0]

	    big_lineages = np.union1d(big_lineages,extant_lineages[sizes[extant_lineages,t+1]>=10])


	    if np.mod(t/2.0,progress_interval)==0 and t/progress_interval/2 <1:
	    	bar.progress(t/progress_interval/2)


	if len(big_lineages)==1:
		big_lineages = np.arange(1,len(lineages))


	return big_lineages,children,sizes,fits

def add_trads(children,sizes,lineage_no):
    if len(children[lineage_no])==0:
        return sizes[lineage_no,:]
    else:
        sumsizess = np.copy(sizes[lineage_no,:])
        for ki in range(len(children[lineage_no])):
            sumsizess += add_trads(children,sizes,children[lineage_no][ki])
        return sumsizess


@st.cache(allow_output_mutation=True)
def plot_trajectories(N,U,s,num_gen,big_lineages,children,sizes):
	chart_data = pd.DataFrame([])
	ki = 1
	arrlist = [list(add_trads(children,sizes,k)) for k in big_lineages[1:]]
	wide_df = pd.DataFrame(arrlist,columns=range(num_gen))
	wide_df['lineage'] = wide_df.index.astype(str)
	chart_data = wide_df.melt(id_vars='lineage',var_name='generation',value_name='number individuals')
	chart_data=chart_data[chart_data['number individuals']>0]


	gen_range = pd.DataFrame({'generation':np.linspace(0,num_gen,100).astype(int)})
	c = alt.Chart(chart_data,width = 400).mark_line().encode(
		x='generation:Q',

		y=alt.Y('number individuals', scale=alt.Scale(domain=[0.0,1.1*N], clamp=True), title='number individuals'),
		color=alt.Color('lineage',legend=None)
	)

	line = alt.Chart(chart_data).mark_line().encode(
		x='generation',

		y=alt.Y('number individuals', scale=alt.Scale(domain=[0.0,1.1*N], clamp=True), title='number individuals'),
		color=alt.Color('lineage',legend=None)
	).interactive()

	nearest = alt.selection(type='single', nearest=True, on='mouseover',fields=['generation'], empty='none')


	slider = alt.binding_range(min=0, max=int(num_gen)-1, step=10, name='generation:')
	selector = alt.selection_single(name="generation", fields=['generation'],
                                bind=slider, init={'generation': 50})

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

	df = pd.DataFrame(sizes[:len(fits),:],columns=range(num_gen),index=fits)
	df['fitness'] = df.index
	df_melt = df.melt(id_vars='fitness',var_name='generation',value_name='number')
	df_melt = df_melt[df_melt.number>0].astype({'generation': 'int32'})

	gen = 50



	bchart = alt.Chart(df_melt,width = 600).mark_bar().encode(
    	x=alt.Y('fitness:Q', bin=alt.Bin(extent=[1, max(df['fitness'])], step=s)),
    	y=alt.Y('number:Q', scale=alt.Scale(domain=(0, N)))).add_selection(nearest).transform_filter(nearest).interactive()


	concat_chart = (layerChart & bchart).configure_axis(labelFontSize=20,titleFontSize=20)
	return concat_chart



if st.sidebar.button('Run simulation'):
	seed = np.random.rand()
	bar = st.progress(0)
	#params = pd.DataFrame({'N':[N],'U':[U],'s':[s],'num_gen':[num_gen],'seed':[seed]})
	#params.to_csv('runparams.csv')


	#params = pd.read_csv('runparams.csv')
	#N,U,s,num_gen,seed = params['N'][0],params['U'][0],params['s'][0],params['num_gen'][0],params['seed'][0]

	big_lineages,children,sizes,fits = runsim(N,U,dfe,num_gen,seed)
	bar.progress(50)
	layerChart = plot_trajectories(N,U,s,num_gen,big_lineages,children,sizes)
	bar.progress(75)
	st.altair_chart(layerChart  )
	bar.progress(100)















