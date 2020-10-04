#!/usr/bin/env python
# encoding: utf-8
"""
EyeLinkSession.py

Created by Tomas Knapen on 2011-04-27.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import os, sys, pickle, math, thread, time
from subprocess import *

import scipy as sp
import scipy.stats as stats
import numpy as np
import matplotlib.pylab as pl
from matplotlib.backends.backend_pdf import PdfPages

from tables import *
import pp

nr_variables = 5

#class for callbacks
class DataContainer:
	def __init__(self, nr_timepoints, nr_simulations, nr_variables):
		self.result_array = np.zeros((nr_simulations, nr_timepoints, nr_variables))
		self.parameter_array = []# [np.zeros((nr_simulations, nr_parameters))]
		self.lock = thread.allocate_lock()
		self.count = 0

	#the callback function
	def save_to_array(self, value):
		# we must use lock here because array stuff is not atomic (?)
		self.lock.acquire()
		self.parameter_array.append( value[0] )
		self.result_array[self.count] = value[1]
		self.lock.release()
		self.count += 1
	
	def save_to_hdf_file(self, hdf5_filename):
		# saving the data
		
		if os.path.isfile(hdf5_filename):
			os.system('rm ' + hdf5_filename)
		h5file = openFile(hdf5_filename, mode = "w", title = " file")
			
		try:
			thisRunGroup = h5file.getNode(where = "/", name = '', classname='Group')
		except NoSuchNodeError:
			# import actual data
			thisRunGroup = h5file.createGroup(where = "/", name = str(sim))
					
		h5file.createArray(thisRunGroup, 'simulation_data', self.result_array, '')
		
		ptd = [(k, np.float64) for k in np.unique(np.concatenate([k.keys() for k in self.parameter_array]))]
		self.parameterTypeDictionary = np.dtype(ptd)
		
		# create a table for the parameters of these runs
		parameterTable = h5file.createTable(thisRunGroup, 'sim_parameters', self.parameterTypeDictionary)
		# fill up the table
		trial = parameterTable.row
		for r in self.parameter_array:
			for par in r.keys():
				trial[par] = r[par]
			trial.append()
		parameterTable.flush()
	
		h5file.close()
		
	def plot_activities(self, plot_file_name, sort_variable = None):
		plot_file = PdfPages(plot_file_name)
		
		if sort_variable == None:
			order = range(len(self.parameter_array))
		else:
			order = np.argsort(np.array([p[sort_variable] for p in self.parameter_array]))
			
			
		for i in order:
			fig = pl.figure(figsize = (15, 6))
			s = fig.add_subplot(211)
			s.set_title('simulation results')
			s.set_xlabel('time [steps]')
			s.set_ylabel('activity strength')
			# s.set_ylim([0,1.3])
			
			pl.plot(self.result_array[i,::10,0], 'r', alpha = 0.75, label = 'H1')
			pl.plot(self.result_array[i,::10,1], 'g', alpha = 0.75, label = 'H2')
			pl.plot(self.result_array[i,::10,4], 'k', alpha = 0.55, label = 'C') #  / np.max(self.result_array[i,::10,4])
			
			leg = s.legend(fancybox = True)
			leg.get_frame().set_alpha(0.5)
			if leg:
				for t in leg.get_texts():
				    t.set_fontsize('small')    # the legend text fontsize
				for l in leg.get_lines():
				    l.set_linewidth(3.5)  # the legend line width
			
			s = fig.add_subplot(212)
			s.set_title('simulation results')
			s.set_xlabel('time [steps]')
			s.set_ylabel('adaptation strength')
			# s.set_ylim([0,1.3])
			
			pl.plot(self.result_array[i,::10,2], 'r--', alpha = 0.25, label = 'A1')
			pl.plot(self.result_array[i,::10,3], 'g--', alpha = 0.25, label = 'A2')
			pl.text(10, 1.0, str(self.parameter_array[i]), fontsize = 8)
			
			leg = s.legend(fancybox = True)
			leg.get_frame().set_alpha(0.5)
			if leg:
				for t in leg.get_texts():
				    t.set_fontsize('small')    # the legend text fontsize
				for l in leg.get_lines():
				    l.set_linewidth(3.5)  # the legend line width
			
			plot_file.savefig()
			pl.close()
		plot_file.close()
		# pl.show()
	

def npS( input, mu ):
	input[input < 0] = 0.0
	input = pow(input,mu['NRa'])/(pow(input,mu['NRa']) + pow(mu['NRs'],mu['NRa']))

def func(t, y, mu):
	import pygsl._numobj
	from pygsl import odeiv, Float
	
	def S( input, NRa, NRs ):
		if input >= 0. :
			 return pow(input,NRa)/(pow(input,NRa) + pow(NRs,NRa))
		else:
			return 0.
	
	dydt = pygsl._numobj.zeros((5,), Float) * 1.0
	
	#defining variables based on indices on y
	H1, H2 = 0,1
	A1, A2 = 2,3
	C = 4
	
	# dydt[H1] = mu['XL'] - (1. + y[A1]) * y[H1] + mu['beta'] * y[A1] - mu['gamma'] * S(y[H2], mu['NRa'], mu['NRs']);
	# dydt[H2] = mu['XR'] - (1. + y[A2]) * y[H2] + mu['beta'] * y[A2] - mu['gamma'] * S(y[H1], mu['NRa'], mu['NRs']);
	
	dydt[H1] = mu['XL'] - (1. + y[A1]) * y[H1] + mu['beta'] * y[A1] - mu['gamma'] * S(y[H2], mu['NRa'], mu['NRs']) - mu['var_inh_infl'] * S(y[C], mu['NRa'], mu['NRs']);
	dydt[H2] = mu['XR'] - (1. + y[A2]) * y[H2] + mu['beta'] * y[A2] - mu['gamma'] * S(y[H1], mu['NRa'], mu['NRs']) - mu['var_inh_infl'] * S(y[C], mu['NRa'], mu['NRs']);
	dydt[A1] = ( -pow(y[A1],mu['exponent']) + ( mu['alpha'] * S(y[H1], mu['NRa'], mu['NRs']) ) ) / mu['tau'];
	dydt[A2] = ( -pow(y[A2],mu['exponent']) + ( mu['alpha'] * S(y[H2], mu['NRa'], mu['NRs']) ) ) / mu['tau'];
	dydt[C] =  (mu['var_inh'] * (S(dydt[H1], mu['NRa'], mu['NRs']) + S(dydt[H2], mu['NRa'], mu['NRs'])) -y[C]) / mu['tau_inh']
	
	return dydt

def run_sim(mu, nr_timepoints, func, npS):
	import pygsl._numobj
	import pygsl
	from pygsl import odeiv, Float
	import numpy
	
	dimension = 5
	step = odeiv.step_rkf45(dimension, func, None, mu)
	control = odeiv.control_y_new(step, 1e-6, 1e-6)
	evolve  = odeiv.evolve(step, control, dimension)
	
	h = 1
	t1 = float(nr_timepoints)
	# initial values - all 0.
	y = pygsl._numobj.array((0.5, 0.5, 0.0, 0.01, 0.0))
	
	op = numpy.zeros((nr_timepoints, dimension))
	iters = 0
	for t in numpy.linspace(0, t1, nr_timepoints):
		t, h, y = evolve.apply(t, t1, h, y)
		op[iters] = y
		y += numpy.concatenate((numpy.random.randn(2) * mu['noise_level'], [0.0, 0.0, 0.0]))
		iters += 1
	
	op = numpy.array(op)
	# naka rushton on activities:
	npS(op[:,0], mu)
	npS(op[:,1], mu)
	# return both output and parameter dictionary
	return [mu, op]

# mu parameters based on dictionary
mu = {'XL': 1.0, 'XR': 1.0, 'beta': 0.24, 'gamma': 3.0, 'exponent': 1.0, 'alpha': 4.0, 'tau': 100.0, 'NRa': 2.0, 'NRs': 1.0, 'noise_level': 0.01, 'var_inh': 120.0, 'tau_inh': 50, 'var_inh_infl': 0.8}

which_var = 'var_inh_infl'
which_values = np.linspace(0.5,1.5,6)

# Create an instance of callback class
nr_simulations = which_values.shape[0]
nr_timepoints = 100000
dc = DataContainer(nr_timepoints = nr_timepoints, nr_simulations = nr_simulations, nr_variables = nr_variables)

# running these in parallel
# Creates jobserver with automatically detected number of workers
job_server = pp.Server(ppservers=())

# Execute the same task with different amount of active workers and measure the time
for index in xrange(nr_simulations):
	mu[which_var] = which_values[index]
	job_server.submit(run_sim, (mu, nr_timepoints, func, npS), callback=dc.save_to_array)
#wait for jobs in all groups to finish 
job_server.wait()

dc.save_to_hdf_file('data/test2.hdf')
# saver.plot_activities('data/test2.pdf', sort_variable = which_var)
