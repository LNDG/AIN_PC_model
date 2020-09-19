#!/usr/bin/env python
# encoding: utf-8
"""
EyeLinkSession.py

Created by Tomas Knapen on 2011-04-27.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import os, sys, pickle, math, thread, time
from subprocess import *
sys.path.append('/home/mpib/kamp/LNDG/Noise_Color_Attractor_Model/data_handling')

import scipy as sp
import scipy.stats as stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as pl
from IPython import embed as shell

from tables import *
import pp
from model_config import *
from DataContainer import DataContainer
from DataAnalyzer import DataAnalyzer

#configurations for the model
nr_variables = 4 #number of variables in the model
# model parameters
nr_timepoints = 60000
simulate=True

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
	
	dydt = pygsl._numobj.zeros((4,), Float) * 1.0
	
	#defining variables based on indices on y
	H1, H2 = 0,1
	A1, A2 = 2,3
	
	dydt[H1] = mu['XL'] - (1. + y[A1]) * y[H1] + mu['beta'] * y[A1] - mu['gamma'] * S(y[H2], mu['NRa'], mu['NRs']) 
	dydt[H2] = mu['XR'] - (1. + y[A2]) * y[H2] + mu['beta'] * y[A2] - mu['gamma'] * S(y[H1], mu['NRa'], mu['NRs'])
	dydt[A1] = ( -pow(y[A1],mu['exponent']) + ( mu['alpha'] * S(y[H1], mu['NRa'], mu['NRs']) ) ) / mu['tau']
	dydt[A2] = ( -pow(y[A2],mu['exponent']) + ( mu['alpha'] * S(y[H2], mu['NRa'], mu['NRs']) ) ) / mu['tau']

	return dydt

def run_sim(mu, nr_timepoints, func, npS):
	import pygsl._numobj
	import pygsl
	from pygsl import odeiv, Float
	import numpy
	
	def S( input, NRa, NRs ):
		if input >= 0. :
			 return pow(input,NRa)/(pow(input,NRa) + pow(NRs,NRa))
		else:
			return 0.00001
	
	dimension = 4
	step = odeiv.step_rkf45(dimension, func, None, mu)
	control = odeiv.control_y_new(step, 1e-6, 1e-6)
	evolve= odeiv.evolve(step, control, dimension)
	
	h = 1
	t1 = float(nr_timepoints)
	# initial values - all 0.
	y = pygsl._numobj.array((0., 0., 0., 0.))
	op = numpy.zeros((nr_timepoints, dimension))
	
	#load noise file 
	noise_dict = {1:'white', 2:'pink', 3:'blue'}
	noise_color = noise_dict[mu['noise_color']]
	noise_file = 'colored_noise/%s_noise.csv' % (noise_color)
	noise = numpy.genfromtxt(noise_file, delimiter=',')
	iters = 0
	for t in numpy.linspace(0, t1, nr_timepoints):
		t, h, y = evolve.apply(t, t1, h, y)
		op[iters] = y
		# add colored noise to instantaneous activity:
		y += numpy.array([noise[0,iters] * mu['noise_level'], noise[1,iters] * mu['noise_level'], 0.0, 0.0])
		iters += 1
	
	op = numpy.array(op)
	
    # naka rushton on activities:
	npS(op[:,0], mu)
	npS(op[:,1], mu)
	# return both output and parameter dictionary
	return [mu, op]

data_dir = 'data/Colored_Noise/'
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
plot_dir = 'plots/Colored_Noise/'
if not os.path.exists(plot_dir):
    os.mkdir(plot_dir)
file_name = data_dir + 'Colored_Noise'
plot_name = plot_dir + 'Colored_Noise'

noise_dict = {1:'white', 2:'pink', 3:'blue'}
corr_res = np.zeros((len(noise_dict), noise_level_range.shape[0]))

for noise_nr, noise_color in noise_dict.items():
	mu['noise_color'] = noise_nr
	which_var = 'noise_level'
	#which_values = inl_range
	rn = noise_color
	
	print 'running simulation with %s noise.' % (noise_color)
	
	# Create an instance of callback class
	nr_simulations = noise_level_range.shape[0]
	dc = DataContainer(file_name + '.hdf5')
	da = DataAnalyzer(dc)
	
	if simulate:
		dc.setup_for_simulation(nr_timepoints = nr_timepoints, nr_simulations = nr_simulations, nr_variables = nr_variables)
		# running these in parallel
		# Creates jobserver with automatically detected number of workers
		job_server = pp.Server(ppservers=())

		# Execute the same task with different amount of active workers and measure the time
		for noise_level in noise_level_range:
			mu['noise_level'] = noise_level
			job_server.submit(run_sim, (mu, nr_timepoints, func, npS), callback=dc.save_to_array)
		#wait for jobs in all groups to finish 
		job_server.wait()
		job_server.destroy()
		
		dc.save_to_hdf_file(run_name = rn)
	
		da.plot_activities(plot_file_name = plot_name + '_act_' + rn + '.pdf', run_name = rn, sort_variable = which_var)
	
	# da.all_time_courses_to_percepts(run_name = rn.replace('.',''), sort_variable = which_var, plot_file_name = file_name + '_' + rn + '.pdf')
	# da.plot_activities(plot_file_name = 'data/act_' + rn + '.pdf', run_name = rn.replace('.',''), sort_variable = which_var)
	# corr_res[i:] = da.correlation_results[:,0]
	
# fig = pl.figure()
# ax = fig.add_subplot(111)	
# cax = ax.imshow(corr_res[i], extent = (which_values[0],which_values[-1],inl_range[0],inl_range[-1]), vmin = 0, vmax = 1)
# cbar = fig.colorbar(cax, ticks=[0, 0.5, 1])
# cbar.ax.set_yticklabels(['0', '0.5', '1'])# vertically oriented colorbar
# ax.set_ylabel('inhibition noise level', fontsize=9)
# ax.set_xlabel('variable inhibition strength', fontsize=9)
# pl.savefig('data/im_' + str(population_noise_level) + '.pdf')

fig = pl.figure()
ax = fig.add_subplot(111)	
for i, cr in enumerate(corr_res):
	pl.plot(inl_range, cr, 'k', alpha = 0.3 + 0.7 * (i / float(corr_res.shape[0])), linewidth = 0.5 + 2.5 * (i / float(corr_res.shape[0])))
pl.plot(inl_range, corr_res.mean(axis = 0), 'k--', linewidth = 5. )
pl.savefig(file_name + '_all.pdf')

from scipy.stats import spearmanr
print spearmanr(inl_range, corr_res.mean(axis = 0))
print [spearmanr(inl_range, c) for c in corr_res]

#pl.show()
shell()