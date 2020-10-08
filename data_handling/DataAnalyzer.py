import os, sys, pickle, math, _thread, time, datetime
import scipy as sp
import scipy.signal as signal
import scipy.stats as stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as pl
from IPython import embed as shell
import seaborn as sns
sns.set()

H1, H2 = 0,1
A1, A2 = 2,3
C = 4

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

class DataAnalyzer(object):
	
	variable_names = ['H1', 'H2', 'A1', 'A2', 'C']
	
	def __init__(self, data_container):
		self.data_container = data_container

	
	def plotNoiseSpectrum(self, Noise, Fsample, axes, noise_freq_range, color="r", xlabel=False):
		"""
		Plots a Single-Sided Amplitude Spectrum of Noise with Sampling Frequency Fsample
		"""
		n = len(Noise) # length of the Noise
		k = np.arange(n)
		T = n/Fsample
		frq = k/T # two sides frequency range
		frq = frq[range(int(n/2))] # one side frequency range

		Y = sp.fft(Noise) # fft computing and normalization
		Y = Y[range(int(n/2))]
		
		axes.plot(frq,np.abs(Y),color, alpha = 0.50) # plotting the spectrum
		
		if xlabel:
			axes.set_xlabel('Frequency (Hz)')
		else:
			axes.tick_params(labelbottom=False)

		axes.set_ylabel('Power', fontsize='small')
		axes.yaxis.set_label_position("right")
		axes.yaxis.tick_right()
	
	def plotHist(self, axes, Noise, color, label, title_suffix):
		axes.set_title('Noise, ' + title_suffix)
		#axes.yaxis.tick_right()
		axes.hist(Noise, bins=50, color=color, label=label, alpha=0.3)
		axes.tick_params(labelleft=False, pad=0.5)
		
	def plot_activities(self, plot_file_name, nr_variables, run_name, noise_freq_range, sort_variable = None):
		self.simulation_parameters, self.simulation_data = self.data_container.data_from_hdf_file(run_name)
		plot_file = PdfPages(plot_file_name)
		
		if sort_variable == None:
			self.order = range(len(self.parameter_array))
		else:
			self.order = np.argsort(np.array([p[sort_variable] for p in self.simulation_parameters]))
			
		for i in self.order:
			fig = pl.figure(figsize = (15, 6))
			nr_grid = int((nr_variables-2)*2)
			grid = pl.GridSpec(nr_grid, 8, wspace=0.4, hspace=0.3)
			nr_datapoints = self.simulation_data.shape[1]
			seconds = nr_datapoints/1000

			# create plot grid
			act_ax = fig.add_subplot(grid[:int(nr_grid/2), :6])
			noise1_ax = fig.add_subplot(grid[int(nr_grid/2):int(nr_grid/2+1),:6])
			spec1_ax = fig.add_subplot(grid[int(nr_grid/2):int(nr_grid/2+1),6:])
			hist = fig.add_subplot(grid[:int(nr_grid/2), 6:])
			noise2_ax = fig.add_subplot(grid[int(nr_grid/2+1):int(nr_grid/2+2),:6])
			spec2_ax = fig.add_subplot(grid[int(nr_grid/2+1):int(nr_grid/2+2),6:])


			if nr_variables == 5:
				noise3_ax = fig.add_subplot(grid[-1,:6])
				spec3_ax = fig.add_subplot(grid[-1,6:])
			
			
			# plot activity
			act_ax.set_title('Simulation Results')
			act_ax.set_ylabel('Activity Strength')
			act_ax.tick_params(labelbottom=False)
			act_ax.plot(self.simulation_data[i,::,0], 'r', alpha = 0.75, label = 'H1')
			act_ax.plot(self.simulation_data[i,::,1], 'g', alpha = 0.75, label = 'H2')
			if nr_variables == 5:
				act_ax2 = act_ax.twinx()
				act_ax2.set_ylabel('Feedback (C)')
				act_ax2.plot(self.simulation_data[i,::,4], 'k', alpha = 0.25, label = 'C')
			type(self.simulation_parameters)
			parameter_string = 'I(1)=' + str(self.simulation_parameters[i][2]) + ' I(2)=' + str(self.simulation_parameters[i][3])
			act_ax.text(0.8, 0.95, parameter_string)

			# nr switches
			nr_switches = self.detectswitches(self.simulation_data[i,:,0], self.simulation_data[i,:,1])
			switch_label = "# Switches = %s" % (nr_switches)
			act_ax.plot([], [], ' ', label=switch_label) 
			
			def setlegend(axes):
				leg = axes.legend(loc='upper right', fancybox = True)
				leg.get_frame().set_alpha(0.5)
				if leg:
					for t in leg.get_texts():
						t.set_fontsize('small')    # the legend text fontsize
					for l in leg.get_lines():
						l.set_linewidth(3.5)  # the legend line width
			setlegend(act_ax)

			# plot noise 
			noise_idx = nr_variables
			noise1_ax.set_ylabel('Noise\nStrength', fontsize='small')
			noise1_ax.tick_params(labelbottom=False)
			noise1_ax.set_ylim(-0.3,0.3)
			noise1_ax.plot(self.simulation_data[i,::,noise_idx], 'r', alpha = 0.50, label = 'Noise H1')
			setlegend(noise1_ax)
			# plot freq spec
			self.plotNoiseSpectrum(self.simulation_data[i,:,noise_idx], Fsample=1000, color='r', axes=spec1_ax, noise_freq_range=noise_freq_range)

			noise2_ax.plot(self.simulation_data[i,::,noise_idx+1], 'g', alpha = 0.50, label = 'Noise H2')
			noise2_ax.set_ylabel('Noise\nStrength', fontsize='small')
			noise2_ax.set_ylim(-0.3,0.3)
			setlegend(noise2_ax)
			
			def setlowerticks(axes):
				axes.set_xlabel('Time [Seconds]')
				axes.set_xticks(np.arange(0,nr_datapoints+1, 1000))
				axes.set_xticklabels(np.arange(0, seconds+1, 1))

			if nr_variables == 5:
				noise2_ax.tick_params(labelbottom=False)
				noise3_ax.plot(self.simulation_data[i,::,noise_idx+2], 'k', alpha = 0.50, label = 'Noise C')
				noise3_ax.set_ylabel('Noise\nStrength', fontsize='small')
				noise3_ax.set_ylim(-0.3,0.3)
				setlegend(noise3_ax)
				setlowerticks(noise3_ax)
				# plot freq spec
				self.plotNoiseSpectrum(self.simulation_data[i,:,noise_idx+1], Fsample=1000, color='g',axes=spec2_ax)
				self.plotNoiseSpectrum(self.simulation_data[i,:,noise_idx+2], Fsample=1000, color='k', axes=spec3_ax, xlabel=True, noise_freq_range=noise_freq_range)
			else:
				setlowerticks(noise2_ax)
				self.plotNoiseSpectrum(self.simulation_data[i,:,noise_idx+1], Fsample=1000, color='g', axes=spec2_ax, xlabel=True, noise_freq_range=noise_freq_range)
			
			# noise plot hist
			title_suffix = 'SD = ' + str(np.round(self.simulation_parameters[i][sort_variable],3))
			self.plotHist(hist, self.simulation_data[i,:,noise_idx], color='r', label='Noise H1', title_suffix=title_suffix)
			self.plotHist(hist, self.simulation_data[i,:,noise_idx+1], color='g', label='Noise H2', title_suffix=title_suffix)
			if nr_variables == 5:
				self.plotHist(hist, self.simulation_data[i,:,noise_idx+2], color='k', label='Noise H2', title_suffix=title_suffix)

			setlegend(hist)

			plot_file.savefig()
			pl.close()
		plot_file.close()
	
	def detectswitches(self, timecourse1, timecourse2):
		diff_tc = timecourse1 - timecourse2
		window = signal.hann(50)
		smooth_dtc = signal.convolve(diff_tc, window, mode='same') / sum(window)
		switch_tc = np.abs(np.diff(np.sign(smooth_dtc)))
		#nr_switches=len(switch_tc[switch_tc!=0])
		nr_switches = np.sum(switch_tc)/2
		return nr_switches
	
	def transition_occurrence_times(self, simulation_data, smoothing_kernel_width = 200):
		#defining variables based on indices on y
		H1, H2 = 0,1
		A1, A2 = 2,3
		
		from scipy import stats, signal
		from scipy.stats import spearmanr
		
		gauss_pdf = stats.norm.pdf( np.linspace(-4, 4, smoothing_kernel_width) )
		gauss_pdf = gauss_pdf / gauss_pdf.sum()
		print(gauss_pdf.shape)

		difference_time_course = simulation_data[:,H1] - simulation_data[:,H2]
		print(difference_time_course.shape)
		smoothed_dtc = np.zeros((difference_time_course.shape[0] + smoothing_kernel_width))
		smoothed_dtc[smoothing_kernel_width/2:-smoothing_kernel_width/2] = signal.fftconvolve(difference_time_course, gauss_pdf, 'same')
		smoothed_dtc = smoothed_dtc[smoothing_kernel_width/2:-smoothing_kernel_width/2]
		#index
		idx = np.array(np.abs((np.diff(np.sign(np.insert(smoothed_dtc,0,0)))/2.0))).astype('bool')
		transition_occurrence_times = np.arange(smoothed_dtc.shape[0])[idx]
		return transition_occurrence_times
	
	def correlate_triad(self, H, C, dur):
		"""docstring for correlate_triad"""
		"TODO eliminate C"
		from scipy.stats import spearmanr
		return [spearmanr(dur, C), spearmanr(H, dur), spearmanr(H, C)] 
	
	def single_time_course_to_percepts(self, simulation_data, smoothing_kernel_width = 200, axScatter = None, axHistx = None, axHisty = None, color_D = 'r', color_H = 'b', value = 0.0, C_sampling_interval = [-25, 50], H_sampling_interval = [50,-25], plot = True, output_raw_data = False, correlate_which_entity = 0):
		"""This method doesn't care at all about transition durations"""
		"TODO eliminate C"
		transition_occurrence_times = self.transition_occurrence_times(simulation_data = simulation_data, smoothing_kernel_width = smoothing_kernel_width)
		transition_occurrence_times = transition_occurrence_times[(transition_occurrence_times > -C_sampling_interval[0]) * (transition_occurrence_times < (simulation_data.shape[0] - C_sampling_interval[1]))]
		
		# separated into on-and off periods:
		transition_occurrence_times_separated = [transition_occurrence_times[::2], transition_occurrence_times[1::2]]
		percept_durations = np.diff(transition_occurrence_times)
		percept_durations_separated = [percept_durations[::2], percept_durations[1::2]]
		
		if transition_occurrence_times_separated[0].shape[0] > 2:
			#C_vals = [[],[]]
			H_vals = [[],[]]
			dur_vals = [[],[]]
			corrs = [[],[]]
			
			for j in [0,1]:
				nr_events = np.min([transition_occurrence_times_separated[j].shape[0], percept_durations_separated[j].shape[0]])
				percept_onset_offsets = np.array([transition_occurrence_times_separated[j][:nr_events],transition_occurrence_times_separated[j][:nr_events] + percept_durations_separated[j][:nr_events]]).T
				# quickly take the mean of the C signal around the time of a transition
				#C_sampling_periods = np.vstack([transition_occurrence_times_separated[j][:nr_events] + C_sampling_interval[0], transition_occurrence_times_separated[j][:nr_events] + C_sampling_interval[1]]).T
				#C_values_per_transition = np.array([simulation_data[c[0]:c[1],4].mean() for c in C_sampling_periods])#[:-1]
				
				H_sampling_periods = np.vstack([percept_onset_offsets[:,0] + H_sampling_interval[0], percept_onset_offsets[:,1] + H_sampling_interval[1]]).T
				H_values_per_percept = np.array([[(simulation_data[h[0]:h[1],0]+simulation_data[h[0]:h[1],1]).mean(), (simulation_data[h[0]:h[1],0]+simulation_data[h[0]:h[1],1]).var(), np.abs((simulation_data[h[0]:h[1],0]-simulation_data[h[0]:h[1],1])).mean(), np.abs((simulation_data[h[0]:h[1],0]-simulation_data[h[0]:h[1],1])).var()] for h in H_sampling_periods])#[:-1]
				
				which_transitions_are_valid = ~np.isnan(C_values_per_transition)
				C_values_per_transition = C_values_per_transition[which_transitions_are_valid]

				percept_durations_for_corr = percept_durations_separated[j][which_transitions_are_valid]
				
				C_vals[j] = C_values_per_transition
				H_vals[j] = H_values_per_percept
				dur_vals[j] = percept_durations_for_corr
				corrs[j] = self.correlate_triad(H_vals[j][:,correlate_which_entity], C_vals[j], dur_vals[j])
				
				if plot:
					if j == 0:
						axScatter.scatter(percept_durations_for_corr, C_values_per_transition, marker = 'o', color = color_D, alpha = 0.1, label = "value = %1.2f, $rho = %1.3f, p = %1.3f" % (value, corrs[j][0][0], corrs[j][0][1]))
					elif j == 1:
						axScatter.scatter(percept_durations_for_corr, C_values_per_transition, marker = '+', color = color_D, alpha = 0.1, label = "value = %1.2f, $rho = %1.3f, p = %1.3f" % (value, corrs[j][0][0], corrs[j][0][1]))
					
					axHistx.hist(percept_durations_for_corr, bins=50, alpha = 0.125, normed = True, histtype = 'stepfilled', linewidth = .5, color = color_D, cumulative = False )
					axHisty.hist(C_values_per_transition[:-1], bins=50, orientation='horizontal', alpha = 0.125, normed = True, histtype = 'stepfilled', linewidth = .5, color = color_D, cumulative = False)
					# axHisty.hist(H_values_per_percept * 1500, bins=200, orientation='horizontal', alpha = 0.125, normed = True, histtype = 'stepfilled', linewidth = 1.5, color = color_H, cumulative = True)

					axHistx.set_xlim( axScatter.get_xlim() )
					axHistx.set_ylim( (-0.1, 1.1) )
					axHisty.set_ylim( axScatter.get_ylim() )
					axHisty.set_xlim( (-0.1, 1.1) )
		else:
			corrs = [[[0, 0], [0, 0], [0, 0]], [[0, 0], [0, 0], [0, 0]]]
			H_values_per_percept, C_values_per_transition, percept_durations = [],[],[]
		
		if output_raw_data:
			return H_vals, C_vals, dur_vals
		else:
			return corrs
	
	def transition_related_averaging_run(self, simulation_data, smoothing_kernel_width = 200, sampling_interval = [-50, 150], plot = True ):
		"""docstring for transition_related_averaging"""
		transition_occurrence_times = self.transition_occurrence_times(simulation_data = simulation_data, smoothing_kernel_width = smoothing_kernel_width)
		# make sure only valid transition_occurrence_times survive
		transition_occurrence_times = transition_occurrence_times[(transition_occurrence_times > -sampling_interval[0]) * (transition_occurrence_times < (simulation_data.shape[0] - sampling_interval[1]))]
		
		# separated into on-and off periods:
		transition_occurrence_times_separated = [transition_occurrence_times[::2], transition_occurrence_times[1::2]]
		
		mean_time_course, std_time_course = np.zeros((2, sampling_interval[1] - sampling_interval[0], 5)), np.zeros((2, sampling_interval[1] - sampling_interval[0], 5))
		if transition_occurrence_times_separated[0].shape[0] > 2:
			
			for k in [0,1]:
				averaging_interval_times = np.array([transition_occurrence_times_separated[k] + sampling_interval[0],transition_occurrence_times_separated[k] + sampling_interval[1]]).T
				interval_data = np.array([simulation_data[avit[0]:avit[1]] for avit in averaging_interval_times])
				mean_time_course[k] = interval_data.mean(axis = 0)
				std_time_course[k] = (interval_data.std(axis = 0) / np.sqrt(interval_data.shape[0]))
			
			if plot:
				f = pl.figure(figsize = (10,8))
				for i in range(simulation_data.shape[1]):
					s = f.add_subplot(simulation_data.shape[1], 1, 1 + i)
					for j in [0,1]:
						pl.plot(np.arange(mean_time_course[j].T[i].shape[0]), mean_time_course[j].T[i], ['r--','b--'][j], linewidth = 2.0 )
						pl.fill_between(np.arange(mean_time_course[j].shape[0]), mean_time_course[j].T[i] + std_time_course[j].T[i], mean_time_course[j].T[i] - std_time_course[j].T[i], ['r','b'][j], alpha = 0.2)
					s.set_title(self.variable_names[i])
				pl.draw()
			
		return (mean_time_course, std_time_course)
	
	def transition_related_averaging(self, run_name = None, sort_variable = None, smoothing_kernel_width = 200, sampling_interval = [-50, 150] ):
		self.ordering(run_name = run_name, sort_variable = sort_variable)
		for d in self.simulation_data:
			self.transition_related_averaging_run(simulation_data = d, smoothing_kernel_width = smoothing_kernel_width, sampling_interval = sampling_interval)
		# pl.show()
	
	def ordering(self, run_name = None, sort_variable = None):
		
		if not hasattr(self, 'simulation_parameters'):
			self.simulation_parameters, self.simulation_data = self.data_container.data_from_hdf_file(run_name)
		
		if sort_variable == None:
			order = range(len(self.simulation_parameters))
		else:
			order = np.argsort(np.array([p[sort_variable] for p in self.simulation_parameters]))
		
		return order
	
	def all_time_courses_to_percepts(self, run_name = None, sort_variable = None, plot_file_name = None):
		# if run_name != None and hasattr(self, 'simulation_parameters'):
		# 	pass
		order = self.ordering(run_name = run_name, sort_variable = sort_variable)
		
		from matplotlib.ticker import NullFormatter
		nullfmt   = NullFormatter()         # no labels
		
		# definitions for the axes
		left, width = 0.1, 0.65
		bottom, height = 0.1, 0.65
		bottom_h = left_h = left+width+0.02
		
		rect_scatter = [left, bottom, width, height]
		rect_histx = [left, bottom_h, width, 0.2]
		rect_histy = [left_h, bottom, 0.2, height]
		
		axScatter = pl.axes(rect_scatter)
		axHistx = pl.axes(rect_histx)
		axHisty = pl.axes(rect_histy)
		
		f1 = pl.figure(1, figsize = (10,10))
		# s = f.add_subplot(111) # , aspect = 'equal')
		# axScatter.set_title('simulation results\ncorrelations between C and percept duration\nfor %s' % sort_variable)
		axScatter.set_xlabel('duration [time - a.u.]')
		axScatter.set_ylabel('C value [a.u.]')
		
		# no labels on histograms
		axHistx.xaxis.set_major_formatter(nullfmt)
		axHisty.yaxis.set_major_formatter(nullfmt)
		
		colors_D = [(1.0-i, i, 0.0) for i in np.linspace(0.0,1.0,len(order))]
		colors_H = [(0.5-i, 0.5-i, i) for i in np.linspace(0.0,0.5,len(order))]
		
		correlation_results = []
		for i in order:
			correlation_results.append(self.single_time_course_to_percepts(self.simulation_data[i], axScatter = axScatter, axHistx = axHistx, axHisty = axHisty, color_D = colors_D[i], color_H = colors_H[i], value = self.simulation_parameters[i][sort_variable], output_raw_data = False))
		correlation_results = np.array(correlation_results)
		
		leg = axScatter.legend(fancybox = True)		
		leg.get_frame().set_alpha(0.5)
		if leg:
			for t in leg.get_texts():
				t.set_fontsize('x-small')    # the legend text fontsize
			for l in leg.get_lines():
				l.set_linewidth(3.5)  # the legend line width
		
		now = datetime.datetime.now()
		pl.savefig(plot_file_name)
		f1.clf()
		pl.close()
		
		f2 = pl.figure(figsize = (8,4))
		s = f2.add_subplot(111) # , aspect = 'equal')
		# s.set_title('simulation results\ncorrelations between C and percept duration\nfor %s' % sort_variable)
		s.set_xlabel('%s' % sort_variable)
		s.set_ylabel('Spearman\'s rho')
		s.plot(self.simulation_parameters[order][sort_variable], correlation_results[:,0,0,0], 'k--', label = 'percept duration / C')
		s.plot(self.simulation_parameters[order][sort_variable], correlation_results[:,0,1,0], 'b--', label = 'percept duration / $\sigma$ H')
		s.plot(self.simulation_parameters[order][sort_variable], correlation_results[:,0,2,0], 'r--', label = '$\sigma$ H / C')
		s.plot(self.simulation_parameters[order][sort_variable], correlation_results[:,1,0,0], 'k', label = 'percept duration / C')
		s.plot(self.simulation_parameters[order][sort_variable], correlation_results[:,1,1,0], 'b', label = 'percept duration / $\sigma$ H')
		s.plot(self.simulation_parameters[order][sort_variable], correlation_results[:,1,2,0], 'r', label = '$\sigma$ H / C')
		s.axhline(0,self.simulation_parameters[order][sort_variable].min()-1,self.simulation_parameters[order][sort_variable].max()+1, linewidth = 0.25)
		leg = s.legend(fancybox = True)
		leg.get_frame().set_alpha(0.5)
		if leg:
			for t in leg.get_texts():
				t.set_fontsize('x-small')    # the legend text fontsize
			for l in leg.get_lines():
				l.set_linewidth(3.5)  # the legend line width
		simpleaxis(s)
		spine_shift(s)
		pl.savefig(os.path.splitext(plot_file_name)[0] + '_corr.pdf')
		f2.clf()
		pl.close()
		self.correlation_results = correlation_results

	def percept_values_for_run(self, run_name = None, sort_variable = None):
		"""docstring for values_for_run"""
		if not hasattr(self, 'simulation_parameters'):
			self.simulation_parameters, self.simulation_data = self.data_container.data_from_hdf_file(run_name)
		order = self.ordering()
		
		return [self.single_time_course_to_percepts(self.simulation_data[i], axScatter = None, axHistx = None, axHisty = None, color_D = 'r', color_H = 'r', value = self.simulation_parameters[i][sort_variable], plot = False, output_raw_data = True) for i in order]
	
	def trial_based_noise_analysis(self, run_name = None, sort_variable = None):
		sim_percept_results = self.percept_values_for_run(run_name = run_name, sort_variable = sort_variable)
		
		H_values = [spr[0] for spr in sim_percept_results]
		C_values = [spr[1] for spr in sim_percept_results]
		percept_values = [spr[2] for spr in sim_percept_results]
		
		
		
		
		
