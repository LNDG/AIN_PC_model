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
from IPython import embed as shell

from tables import *
import pp

from DataContainer import DataContainer
from DataAnalyzer import DataAnalyzer

def simpleaxis(ax):
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	
def spine_shift(ax, shift = 10):
	for loc, spine in ax.spines.iteritems():
		if loc in ['left','bottom']:
			spine.set_position(('outward', shift)) # outward by 10 points
		elif loc in ['right','top']:
			spine.set_color('none') # don't draw spine
		else:
			raise ValueError('unknown spine location: %s'%loc)


nr_variables = 5

def npS( input, mu ):
	input[input < 0] = 0.0
	input = pow(input,mu['NRa'])/(pow(input,mu['NRa']) + pow(mu['NRs'],mu['NRa']))

def func(t, y, mu):
	import pygsl._numobj
	from pygsl import odeiv, Float
	from numpy.random import randn
	
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
	
	dydt[H1] = mu['XL'] - (1. + y[A1]) * y[H1] + mu['beta'] * y[A1] - mu['gamma'] * (mu['var_inh_noise_infl'] * y[C] + 1.0) * S(y[H2], mu['NRa'], mu['NRs']);# - mu['var_inh_infl'] * S(y[C], mu['NRa_var_inh'], mu['NRs_var_inh']);
	dydt[H2] = mu['XR'] - (1. + y[A2]) * y[H2] + mu['beta'] * y[A2] - mu['gamma'] * (mu['var_inh_noise_infl'] * y[C] + 1.0) * S(y[H1], mu['NRa'], mu['NRs']);# - mu['var_inh_infl'] * S(y[C], mu['NRa_var_inh'], mu['NRs_var_inh']);
	dydt[A1] = ( -pow(y[A1],mu['exponent']) + ( mu['alpha'] * S(y[H1], mu['NRa'], mu['NRs']) ) ) / mu['tau'];
	dydt[A2] = ( -pow(y[A2],mu['exponent']) + ( mu['alpha'] * S(y[H2], mu['NRa'], mu['NRs']) ) ) / mu['tau'];
	dydt[C] =  (randn(1) * mu['var_inh_noise_level'] - y[C]) / mu['tau_inh']
	
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
		# add noise to instantaneous activity:
		# y += numpy.concatenate((numpy.random.randn(2) * mu['noise_level'], [0.0, 0.0, 0.0]))
		# add noise to novel interaction
		# y += numpy.concatenate(([0.0, 0.0, 0.0, 0.0], numpy.random.randn(1) * mu['noise_level']))
		# add noise to activities and to novel interaction
		# y += numpy.array([numpy.random.randn(1) * mu['noise_level'] * mu['var_inh_noise_infl']/y[4], numpy.random.randn(1) * mu['noise_level'] * mu['var_inh_noise_infl']/y[4], 0.0, 0.0, y[4] * numpy.random.randn(1) * mu['var_inh_noise_level']])
		# add noise only to novel interaction, but graded by the inverse of its value.
		# y += numpy.concatenate(([0.0, 0.0, 0.0, 0.0], numpy.random.randn(1) * mu['noise_level']))
		# add noise to both populations and transient signal
		y += numpy.array([numpy.random.randn(1) * mu['noise_level'], numpy.random.randn(1) * mu['noise_level'], 0.0, 0.0, 0.0])
		
		iters += 1
	
	op = numpy.array(op)
	# naka rushton on activities:
	npS(op[:,0], mu)
	npS(op[:,1], mu)
	# return both output and parameter dictionary
	return [mu, op]

# mu parameters based on dictionary
mu = {'XL': 1.0, 'XR': 1.0, 'beta': 0.24, 'gamma': 3.0, 'exponent': 1.0, 'alpha': 3.0, 'tau': 100.0, 'NRa': 2.0, 'NRs': 1.0, 'noise_level': 0.0025, 'var_inh': 120.0, 'tau_inh': 50, 'var_inh_infl': 0.8, 'NRa_var_inh': 3.0, 'NRs_var_inh': 1.0, 'var_inh_noise_level': 0.005, 'var_inh_noise_infl': 0.0}
nr_timepoints = 20000
file_name = 'data/C_noise_no_transition'

corr_res = np.zeros((4,4,8,6))
pnl_range = np.linspace(0.0001, 0.0003, corr_res.shape[0])
inl_range = np.linspace(0.001, 0.01, corr_res.shape[1])

simulate = True 

for i, population_noise_level in enumerate(pnl_range):
	for j, inhibition_noise_level in enumerate(inl_range):
		mu['var_inh_noise_level'] = inhibition_noise_level
		mu['noise_level'] = population_noise_level
		
		rn = 'pnl' + '_' + str(population_noise_level) + '_inl' + '_' + str(inhibition_noise_level)
		
		which_var = 'var_inh_infl'
		which_values = np.linspace(0.00,7.5,corr_res.shape[2])
		
		# Create an instance of callback class
		nr_simulations = which_values.shape[0]
		dc = DataContainer(file_name + '.hdf5')  
		da = DataAnalyzer(dc)
		if simulate:
			dc.setup_for_simulation(nr_timepoints = nr_timepoints, nr_simulations = nr_simulations, nr_variables = nr_variables)
			# running these in parallel
			# Creates jobserver with automatically detected number of workers
			job_server = pp.Server(ppservers=())
		
			# Execute the same task with different amount of active workers and measure the time
			for index in xrange(nr_simulations):
				mu[which_var] = which_values[index]
				job_server.submit(run_sim, (mu, nr_timepoints, func, npS), callback=dc.save_to_array)
			#wait for jobs in all groups to finish 
			job_server.wait()
			job_server.destroy()
		
			dc.save_to_hdf_file(run_name = rn.replace('.',''))
			da.plot_activities(plot_file_name = file_name + '_act_' + rn + '.pdf', run_name = rn.replace('.',''), sort_variable = which_var)
		
		da.all_time_courses_to_percepts(run_name = rn.replace('.',''), sort_variable = which_var, plot_file_name = file_name + '_' + rn + '.pdf')
		corr_res[i,j,:,:] = da.correlation_results
		
	# fig = pl.figure()
	# ax = fig.add_subplot(111)	
	# cax = ax.imshow(corr_res[i], extent = (which_values[0],which_values[-1],inl_range[0],inl_range[-1]), vmin = 0, vmax = 1)
	# cbar = fig.colorbar(cax, ticks=[0, 0.5, 1])
	# cbar.ax.set_yticklabels(['0', '0.5', '1'])# vertically oriented colorbar
	# ax.set_ylabel('inhibition noise level', fontsize=9)
	# ax.set_xlabel('variable inhibition strength', fontsize=9)
	# pl.savefig('data/im_' + str(population_noise_level) + '.pdf')
	# pl.close()

# # for a run of 1x10:
# cr_m = corr_res.squeeze().mean(axis = 0)
# cr_s = corr_res.squeeze().std(axis = 0) / math.sqrt(10)
# 
# 
# f2 = pl.figure(figsize = (8,4))
# s = f2.add_subplot(111) # , aspect = 'equal')
# # s.set_title('simulation results\ncorrelations between C and percept duration\nfor %s' % sort_variable)
# s.set_xlabel('Strength of transient signal [C] influence')
# s.set_ylabel('Spearman\'s $\rho$')
# s.axhline(0,-0.5,8.0, linewidth = 0.25)
# s.plot(which_values, cr_m[:,0], 'k--', label = 'percept duration / C')
# s.plot(which_values, cr_m[:,2], 'b--', label = 'percept duration / $\sigma$ H')
# s.plot(which_values, cr_m[:,4], 'r--', label = '$\sigma$ H / C')
# 
# pl.fill_between(which_values, cr_m[:,0] + cr_s[:,0], cr_m[:,0] - cr_s[:,0], color = 'k', alpha = 0.2)
# pl.fill_between(which_values, cr_m[:,2] + cr_s[:,2], cr_m[:,2] - cr_s[:,2], color = 'b', alpha = 0.2)
# pl.fill_between(which_values, cr_m[:,4] + cr_s[:,4], cr_m[:,4] - cr_s[:,4], color = 'r', alpha = 0.2)
# s.axis([-0.5,8.0, -1, 1])
# 
# leg = s.legend(fancybox = True)
# leg.get_frame().set_alpha(0.5)
# if leg:
# 	for t in leg.get_texts():
# 	    t.set_fontsize('x-small')    # the legend text fontsize
# 	for l in leg.get_lines():
# 	    l.set_linewidth(3.5)  # the legend line width
# simpleaxis(s)
# spine_shift(s)
# pl.savefig(file_name + '_corr.pdf')
# 

# shell()