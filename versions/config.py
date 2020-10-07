# chose model 
model_name = 'gaba'

# Define default parameter
default_params = {'nr_timepoints':3000, 'XL': 1.0, 'XR': 1.0, 'gamma': 3.25, 'exponent': 1.0, 'alpha': 1.75, 'tau': 100.0, 'NRa': 2.0, 'NRs': 1.0, 
                'noise_lowcut': 5, 'noise_highcut': 300, 'nr_color': 1, 'noise_level':0.001}

def create_params(**kwargs):
	"""
	Function to create model params
	""" 
	params = default_params
	for key, item in kwargs.items():
		params.update({key: item})
	return params

# set up over which variables to iterate
import numpy as np
param_grid = {'XL':np.round(np.arange(0.4,1.4,0.2),2)}

# set parallel variable 
parallel_var = 'noise_level'
parallel_range=np.arange(0.0025, 0.1, 0.05) # defines the range of noise level


