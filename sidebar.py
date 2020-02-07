import numpy as np
import streamlit as st
import pandas as pd 
import altair as alt 
import scipy.special

def run_sidebar():

	N = 10**st.sidebar.slider('choose a population size',min_value=2.0,max_value=4.0,step=0.01,value=3.2,format='10E%f')

	beneficial = st.sidebar.checkbox('beneficial mutations',value=True)
	deleterious = st.sidebar.checkbox('deleterious mutations',value=False)

	if beneficial:
		Ub = 10**st.sidebar.slider('choose a beneficial mutation rate',min_value=-5.0,max_value=-2.0,step=0.01,value=-3.0,format='10E%f')
	else:
		Ub = 0
	if deleterious:
		Ud = 10**st.sidebar.slider('choose a deleterious mutation rate',min_value=-5.0,max_value=-2.0,step=0.01,value=-3.0,format='10E%f')
	else:
		Ud = 0

	if beneficial:
		sb = 10**st.sidebar.slider('choose an average beneficial effect size',min_value=-4.0,max_value=-1.0,step=0.01,value=-2.0,format='10E%f')
	else:
		sb = 0

	if deleterious:
		sd = 10**st.sidebar.slider('choose an average deleterious effect size',min_value=-4.0,max_value=-1.0,step=0.01,value=-2.0,format='10E%f')
	else:
		sd = 0

	if beneficial:
		dfeb_option = st.sidebar.selectbox('choose a distribution of fitness effects of beneficial mutations',pd.Series(['single effect','exponential','gamma']))

		if dfeb_option == 'gamma':

			gamma_shape_b = st.sidebar.slider('choose a (beneficial) gamma shape parameter',min_value=1.0,max_value=32.0,step=1.0,value=1.0)
	else:
		dfeb_option = 'none'

	if deleterious:
		dfed_option = st.sidebar.selectbox('choose a distribution of fitness effects of deleterious mutations',pd.Series(['single effect','exponential','gamma']))

		if dfed_option == 'gamma':

			gamma_shape_d = st.sidebar.slider('choose a (deleterious) gamma shape parameter',min_value=1.0,max_value=32.0,step=1.0,value=1.0)
	else:
		dfed_option = 'none'


	num_gen = st.sidebar.slider('how many generations?',min_value=20,max_value=10000,step=10,value=1000)
	assay_interval = st.sidebar.slider('choose the assay interval',min_value=1,max_value=int(num_gen/10),step=1,value=10)


	sbvals = np.linspace(0,4*sb,100)
	sdvals = -np.linspace(4*sd,0,100)


	if beneficial:
		if dfeb_option == 'exponential':
			dfe_source = pd.DataFrame({
			  'fitness effect': sbvals,
			  'probability': Ub/(Ub+Ud)/sb*np.exp(-sbvals/sb)
			})
			dfe_b = lambda : np.random.exponential(sb)


			cb = alt.Chart(dfe_source,width=300).mark_area().encode(
		    x='fitness effect',
		    y='probability'

			)

		elif dfeb_option == 'gamma':
			beta_b = gamma_shape_b/sb
			dfe_source = pd.DataFrame({
			  'fitness effect': sbvals,
			  'probability': Ub/(Ub+Ud)*beta_b**gamma_shape_b/scipy.special.gamma(gamma_shape_b)*sbvals**(gamma_shape_b-1)*np.exp(-beta_b*sbvals)
			})

			dfe_b = lambda : np.random.gamma(gamma_shape_b,1/beta_b)

			cb = alt.Chart(dfe_source,width=300).mark_area().encode(
		    x='fitness effect',
		    y='probability'

			)
			


		elif dfeb_option == 'single effect':
			delta = 0.02
			dfe_source = pd.DataFrame({
			  'fitness effect': [sb*(1-delta),sb*(1+delta)],
			  'probability': [Ub/(Ub+Ud)/2.0/delta/sb,Ub/(Ub+Ud)/2.0/delta/sb]
			})
			dfe_b = lambda : sb


			cb = alt.Chart(dfe_source,width=300).mark_area().encode(
		    	x=alt.X('fitness effect', scale=alt.Scale(domain=[-4*sd,4*sb])),
		    	y='probability'
			)
	else:
		dfe_b = lambda : 0


	if deleterious:
		if dfed_option == 'exponential':
			dfed_source = pd.DataFrame({
			  'fitness effect': sdvals,
			  'probability': Ud/(Ub+Ud)/sd*np.exp(sdvals/sd)
			})
			dfe_d = lambda : np.random.exponential(sd)


			cd = alt.Chart(dfed_source,width=300).mark_area(color='red').encode(
		    	x='fitness effect',
		    	y='probability'

			)

		elif dfed_option == 'gamma':
			beta_d = gamma_shape_d/sd
			dfed_source = pd.DataFrame({
			  'fitness effect': sdvals,
			  'probability': Ud/(Ub+Ud)*beta_d**gamma_shape_d/scipy.special.gamma(gamma_shape_d)*(-sdvals)**(gamma_shape_d-1)*np.exp(beta_d*sdvals)
			})
			dfe_d = lambda : np.random.gamma(gamma_shape_d,1/beta_d)

			cd = alt.Chart(dfed_source,width=300).mark_area(color='red').encode(
		    x='fitness effect',
		    y='probability'

			)
			
		elif dfed_option == 'single effect':
			delta = 0.02
			dfe_source = pd.DataFrame({
			  'fitness effect': [-sd*(1-delta),-sd*(1+delta)],
			  'probability': [Ud/(Ub+Ud)/2.0/delta/sd,Ud/(Ub+Ud)/2.0/delta/sd]
			})
			dfe_d = lambda : sd

			cd = alt.Chart(dfe_source,width=300).mark_area(color='red').encode(
		    	x=alt.X('fitness effect', scale=alt.Scale(domain=[-4*sd,4*sb])),
		    	y='probability'
			)
	else:
		dfe_d = lambda : 0


	if beneficial and deleterious:
		st.sidebar.altair_chart(cb+cd)
	elif beneficial:
		st.sidebar.altair_chart(cb)
	elif deleterious:
		st.sidebar.altair_chart(cd)

	
	return N,Ub,Ud,sb,sd,dfe_b,dfe_d,num_gen,assay_interval

