import os, sys, pickle, math, threading, time
from subprocess import *
import numpy as np
from tables import *
import pygsl._numobj
import pygsl
from pygsl import odeiv, Float
from joblib import Parallel, delayed
from data_handling.DataContainer import DataContainer
from colored_noise.noise import *

def npS(input, params):
	"""""
	Naka Rushton Function
	"""""
	input[input < 0] = 0.0
	input = pow(input,params['NRa'])/(pow(input,params['NRa']) + pow(params['NRs'],params['NRa']))

def integrate_model(func, params):
	""""
	Function to integrate the model ode using the pygsl library
	""""
	# set up pygsl
	step = odeiv.step_rkf45(params['dimension'], func, None, params) # Embedded 4th order Runge-Kutta-Fehlberg method with 5th order error estimate. 
	control = odeiv.control_y_new(step, 1e-6, 1e-6)
	evolve= odeiv.evolve(step, control, params['dimension'])
	
	h = 1
	t1 = float(params['nr_timepoints'])

	y = init_values()
	op = np.zeros((params['nr_timepoints'], params['dimension']))
	
	# init and load noise traces
	noise = noise()
	noise_tc = np.array([[]])
	iters = 0
	for t in np.linspace(0, t1, params['nr_timepoints']):
		t, h, y = evolve.apply(t, t1, h, y)
		op[iters] = y
		# add colored noise to instantaneous activity:
		noise_step = noise.create_step(iters)
		y += noise_step
		if iters == 0:
			noise_tc = noise_step[np.newaxis,:params['nr_noise_tc']]
		else:
			noise_tc = np.concatenate((noise_tc, noise_step[np.newaxis,[0,1]]), axis=0)		
		iters += 1	
	op = np.array(op)	
	
	# naka rushton on activities:
	for act_tc in params['nr_act_tc']:
		npS(op[:,act_tc], params)
	# join noise values to output 
	op = np.concatenate((op, noise_tc), axis=1)
	return [params, op]

def params_generator(params, variable, variable_range):
	""""
	Generator to update model parameters in each iteration
	""""
	assert (variable in params.keys()), 'Variable not found.'	
	for v in variable_range:
		yield params.update({variable: v})		

def run_parallel_integration(func, params, variable, variable_range, hdf5file, hdf5node):
	"""""
	Function to run simulation in parallel over a range of values of one variable 
	"""""
	nr_simulations = variable_range.shape[0]
	# Create an instance of data container class
	dc = DataContainer(hdf5file + '.hdf5')
	dc.setup_for_simulation(nr_timepoints = params['nr_timepoints'], nr_simulations = nr_simulations, 
											nr_variables = params['dimension'], nr_noise = params['nr_noise_tc'])
	# running these in parallel
	Parallel(n_jobs=nr_simulations)(delayed(dc.save_to_array)(integrate_model(func, params)) for _ in params_generator(params, variable, variable_range))
	dc.save_to_hdf_file(run_name = hdf5node)
	