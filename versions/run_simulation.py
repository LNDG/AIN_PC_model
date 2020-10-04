# script to run simulation
from versions.models import *
from versions.parellel_integration import *

simulate = True # set to false to only generate new plots

# define directories for simulated data and plots
data_dir = 'data/testing'
if not os.path.exists(data_dir):
	os.mkdir(data_dir)
plot_dir = 'plots/testing'
if not os.path.exists(plot_dir):
	os.mkdir(plot_dir)

# set file_name of hdf5 data file
hdf5file = data_dir + 'Colored_Noise_Trial_' + noise_color + '_' + str(noise_lowcut) + '-' + str(noise_highcut)
if simulate and os.path.exists(hdf5file+'.hdf5'):
	print(f'Deleting {hdf5file}.hdf5'))
	os.remove(hdf5file+'.hdf5')

# set plot file name
pdffile = plot_dir + 'Colored_Noise_Trial_' + noise_color + '_' + str(noise_lowcut) + '-' + str(noise_highcut)

# set up model
model = gaba()

# set up for loop over run variable
run_range = np.arange(0.8, 1.3, 0.1)
run_variable = 'XL'

def create_runname(**kwargs): 
	for 

for run_value in run_range:
	model.params[run_variable] = run_value
	run_name = f'{run_variable}_{str(run_value).replace('.','')}'
	hdfnode = run_name

	print(f'running simulation for run variable {run_variable} = {run_value}')
	
	# variable for parallel simulation
	parallel_var = 'noise_level'
	parallel_range=np.arange(0.0025, 0.1, 0.0025) # defines the range of noise level

	# Create an instance of callback class
	nr_simulations = noise_level_range.shape[0]
	dc = DataContainer(hdf5file + '.hdf5')

	# Run simulation
	if simulate:
		run_parallel_integration(model, parallel_var, parallel_range, hdf5file, hdf5node)
	
	# Create instance of data analyzer
	da = DataAnalyzer(dc)
	noise_freq_range = [model.noise_params['lowcut'], model.noise_params['highcut']]
	da.plot_activities(plot_file_name = pdffile + '_' + run_name + '.pdf', nr_variables = model.params['dimension'], run_name = run_name, sort_variable = parallel_var, noise_freq_range=noise_freq_range)
