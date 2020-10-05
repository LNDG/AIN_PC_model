"""
EyeLinkSession.py

Created by Tomas Knapen on 2011-04-27.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import os, sys, pickle, math, threading, time
from subprocess import *
sys.path.append('/home/mpib/kamp/LNDG/Noise_Color_Attractor_Model/data_handling')

import numpy as np
from tables import *

def npS(input, params):
	"""""
	Naka Rushton Function
	"""""
	input[input < 0] = 0.0
	input = pow(input,params['NRa'])/(pow(input,params['NRa']) + pow(params['NRs'],params['NRa']))

def integrate_model(model):
	import pygsl._numobj
	import pygsl
	from pygsl import odeiv, Float
	import numpy as np

	step = odeiv.step_rkf45(model.params['dimension'], model.func, None) # Embedded 4th order Runge-Kutta-Fehlberg method with 5th order error estimate. 
	control = odeiv.control_y_new(step, 1e-6, 1e-6)
	evolve= odeiv.evolve(step, control, model.params['dimension'])
	
	h = 1
	t1 = float(model.params['nr_timepoints'])

	y = model.init_values
	op = np.zeros((model.params['nr_timepoints'], model.noise_params['dimension']))
	
	noise = model.noise
	noise_tc = np.array([[]])
	iters = 0
	for t in np.linspace(0, t1, model.params['nr_timepoints']):
		t, h, y = evolve.apply(t, t1, h, y)
		op[iters] = y
		# add colored noise to instantaneous activity:
		y += model.create_noise_step(iters)
		if iters == 0:
			noise_tc = noise_step[np.newaxis,[0,1]]
		else:
			noise_tc = np.concatenate((noise_tc, noise_step[np.newaxis,[0,1]]), axis=0)		
		iters += 1
	
	op = np.array(op)
	
	# naka rushton on activities:
	npS(op[:,0], model.params)
	npS(op[:,1], model.params)
	# join noise values to output 
	op = np.concatenate((op, noise_tc), axis=1)
	# return both output and parameter dictionary
	return [model.params, op]

def params_generator(model, variable, variable_range):
	if variable in model.params.keys():
		params = model.params
	elif variable in model.noise_params.keys():
		params = model.noise_params
	else: raise Exception('Variable not found.')
	
	for v in variable_range:
		params.update({variable: v})
		yield 

def run_parallel_integration(model, variable, variable_range, hdf5file, hdf5node):
	"""""
	Function to run simulation in parallel over a range of values of one variable 
	"""""
	from joblib import Parallel, delayed
	from DataContainer import DataContainer

	nr_simulations = variable_range.shape[0]
	# Create an instance of data container class
	dc = DataContainer(hdf5file + '.hdf5')
	dc.setup_for_simulation(nr_timepoints = model.params['nr_timepoints'], nr_simulations = nr_simulations, 
											nr_variables = model.params['dimension'], nr_noise = model.noise_params['nr_noise_tc'])
	# running these in parallel
	Parallel(n_jobs=nr_simulations)(delayed(dc.save_to_array)(integrate_model(model)) for _ in params_generator(model, variable, variable_range))
	dc.save_to_hdf_file(run_name = hdf5node)
	