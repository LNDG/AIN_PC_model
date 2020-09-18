#!/usr/bin/env python
# encoding: utf-8
"""
EyeLinkSession.py

Created by Tomas Knapen on 2011-04-27.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import os, sys, pickle, math, thread, time, datetime
from subprocess import *

import scipy as sp
import scipy.stats as stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as pl

# from IPython import embed as shell

from tables import *

#class for callbacks
class DataContainer(object):
	def __init__(self, hdf5_filename):
		self.hdf5_filename = hdf5_filename
		
	def setup_for_simulation(self, nr_timepoints, nr_simulations, nr_variables):
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
	
	def save_to_hdf_file(self, run_name, add_to_file = True):
		# saving the data
		
		if os.path.isfile(self.hdf5_filename) and add_to_file:
			# os.system('rm ' + self.hdf5_filename)
			h5file = open_file(self.hdf5_filename, mode = "r+", title = "simulation results file")
		elif os.path.isfile(self.hdf5_filename) and not add_to_file:
			os.system('rm ' + self.hdf5_filename)
			h5file = open_file(self.hdf5_filename, mode = "w", title = "simulation results file")
		else:
			h5file = open_file(self.hdf5_filename, mode = "w", title = "simulation results file")
			
		try:
			thisRunGroup = h5file.get_node(where = "/", name = run_name, classname='Group')
		except NoSuchNodeError:
			# import actual data
			now = datetime.datetime.now()
			thisRunGroup = h5file.create_group("/", run_name, run_name + ' created at ' + now.strftime("%Y-%m-%d_%H.%M.%S"))
		h5file.create_array(thisRunGroup, 'simulation_data', self.result_array, '')
		
		ptd = [(k, np.float64) for k in np.unique(np.concatenate([k.keys() for k in self.parameter_array]))]
		self.parameterTypeDictionary = np.dtype(ptd)
		
		# create a table for the parameters of these runs
		parameterTable = h5file.create_table(thisRunGroup, 'simulation_parameters', self.parameterTypeDictionary)
		# fill up the table
		trial = parameterTable.row
		for r in self.parameter_array:
			for par in r.keys():
				trial[par] = r[par]
			trial.append()
		parameterTable.flush()
	
		h5file.close()
	
	def data_from_hdf_file(self, run_name):
		if not os.path.isfile(self.hdf5_filename):
			print self.hdf5_filename + ' is  not a file'
		self.h5file = open_file(self.hdf5_filename, mode = "r")
		
		try:
			thisRunGroup = self.h5file.get_node(where = "/", name = run_name, classname='Group')
		except NoSuchNodeError:
			# import actual data
			print run_name + ' is not a run in ' + self.hdf5_filename
			return (None, None)
		
		simulation_parameters, simulation_data = thisRunGroup.simulation_parameters.read(), thisRunGroup.simulation_data.read()
		# shell()
		self.h5file.close()
		
		return (simulation_parameters, simulation_data)
		
