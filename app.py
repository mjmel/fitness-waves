import streamlit as st
# To make things easier later, we're also importing numpy and pandas for
# working with sample data.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import time 
import runsim as rs
import altair as alt
import scipy.special
#%%



#%%

st.title('interactive evolution simulations')
st.text('Matthew Melissa, PhD Candidate at Harvard University')

st.markdown('evolution involves a complex interplay of **mutation**, **selection** and **genetic drift**') 

popsize = 10**st.sidebar.slider('choose a population size',min_value=2.0,max_value=4.0,step=0.01,value=3.0,format='10E%f')
mutation_rate = 10**st.sidebar.slider('choose a mutation rate',min_value=-5.0,max_value=-2.0,step=0.01,value=-3.0,format='10E%f')


dfe_option = st.sidebar.selectbox('choose a distribution of fitness effects of mutations',pd.Series(['exponential','gamma']))

if dfe_option == 'gamma':

	gamma_shape = st.sidebar.slider('choose a gamma shape parameter',min_value=1.0,max_value=32.0,step=1.0,value=1.0)



effect_size = 10**st.sidebar.slider('choose an average effect size',min_value=-4.0,max_value=-1.0,step=0.01,value=-2.0,format='10E%f')

svals = np.linspace(0,4*effect_size,100)

if dfe_option == 'exponential':
	dfe_source = pd.DataFrame({
	  'fitness effect': svals,
	  'probability': 1/effect_size*np.exp(-svals/effect_size)
	})
elif dfe_option == 'gamma':
	beta = gamma_shape/effect_size
	dfe_source = pd.DataFrame({
	  'fitness effect': svals,
	  'probability': beta**gamma_shape/scipy.special.gamma(gamma_shape)*svals**(gamma_shape-1)*np.exp(-beta*svals)
	})

#st.sidebar.area_chart(dfe_source).encode(x='fitness effect',y='probability')
c = alt.Chart(dfe_source,width=300).mark_area().encode(
    x='fitness effect',
    y='probability'

)
st.sidebar.altair_chart(c)

def add_trads(lineage_no):
    if len(children[lineage_no])==0:
        return sizes[lineage_no,:]
    else:
        sumsizess = np.copy(sizes[lineage_no,:])
        for ki in range(len(children[lineage_no])):
            sumsizess += add_trads(children[lineage_no][ki])
        return sumsizess


if st.sidebar.button('Run simulation'):
	num_gen = 500




	with st.spinner('Evolving...'):
		big_lineages,children,sizes = rs.runsim(popsize,mutation_rate,effect_size,num_gen)


	chart_data = pd.DataFrame([])
	for k in big_lineages[1:]:
	    chart_data = chart_data.append(pd.Series(add_trads(k)),ignore_index=True)

	aa = st.line_chart(chart_data.T)
	st.text('number of generations')

	st.markdown('each colored line represents the number of individuals with a certain **mutation** in the population')

	st.success('Done!')
	##plt.figure()
	##plt.plot(chart_data.T)

	st.pyplot()

chart_data = pd.DataFrame(
     np.random.randn(200, 3),
     )



x = np.arange(100)
source = pd.DataFrame({
  'x': x,
  'f(x)': np.sin(x / 5)
})

c = alt.Chart(source).mark_line().encode(
    x='x',
    y='f(x)'
)
#%st.altair_chart(c)
