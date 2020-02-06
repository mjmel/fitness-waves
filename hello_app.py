import streamlit as st
import numpy as np

@st.cache(hash_funcs={np.ufunc:id})
def get_num(a):
	return np.exp(a)
st.text(get_num(2))

st.title('Hello world')
