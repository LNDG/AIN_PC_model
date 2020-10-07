import pygsl._numobj
from pygsl import odeiv, Float
import numpy as np
import config

def S( input, NRa, NRs ):
	if input >= 0. :
		return pow(input,NRa)/(pow(input,NRa) + pow(NRs,NRa))
	else:
		return 0.

def set_up(**grid): 
	params = config.create_params(dimension=4, nr_noise_tc=2, nr_act_tc=2)
	init_values = pygsl._numobj.array((0.5,0.5,0,0))
	return params, init_values

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
