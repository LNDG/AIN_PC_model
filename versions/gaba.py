import pygsl._numobj
from pygsl import odeiv, Float
import numpy as np

def default_params():
	params = {'nr_timepoints':3000, 'XL': 1.0, 'XR': 1.0, 'gamma': 3.25, 'exponent': 1.0, 'alpha': 1.75, 'tau': 100.0, 'NRa': 2.0, 'NRs': 1.0, 'noise_lowcut': 5, 'noise_highcut': 300, 'nr_color': 1, 'noise_level':0.001}
	return params

def create_params(**kwargs): 
	params = default_params()
	params.update({'dimension':4, 'nr_noise_tc':2, 'nr_act_tc':2})
	for key, item in kwargs.items():
		params.update({key: item})
	return params

def S( input, NRa, NRs ):
	if input >= 0. :
		return pow(input,NRa)/(pow(input,NRa) + pow(NRs,NRa))
	else:
		return 0.

def func(t, y, params):
	dydt = pygsl._numobj.zeros((4,), Float) * 1.0
	
	#defining variables based on indices on y
	H1, H2 = 0,1
	A1, A2 = 2,3
	
	dydt[H1] = params['XL'] - (1. + y[A1]) * y[H1] - params['gamma'] * S(y[H2], params['NRa'], params['NRs']) 
	dydt[H2] = params['XR'] - (1. + y[A2]) * y[H2] - params['gamma'] * S(y[H1], params['NRa'], params['NRs'])
	dydt[A1] = ( -pow(y[A1],params['exponent']) + ( params['alpha'] * S(y[H1], params['NRa'], params['NRs']) ) ) / params['tau']
	dydt[A2] = ( -pow(y[A2],params['exponent']) + ( params['alpha'] * S(y[H2], params['NRa'], params['NRs']) ) ) / params['tau']

	return dydt 

def init_values():
	# initial values
	init_values = pygsl._numobj.array((0.5, 0.5, 0., 0.))
	return init_values
	
