# script to run simulation
from versions.models import *
from versions.parellel_integration import *
# choose model 
from versions.gaba import *
from data_handling.DataContainer import DataContainer
from data_handling.DataAnalyzer import DataAnalyzer

simulate = True # set to false to only generate new plots

# define directories for simulated data and plots
data_dir = os.path.join('data','testing')
if not os.path.exists(data_dir):
	os.makedirs(data_dir)
plot_dir = os.path.join('plots','testing')
if not os.path.exists(plot_dir):
	os.makedirs(plot_dir)

# set up model parameter
color_dict = {1:'white', 2:'pink', 3:'blue'}
noise_color = color_dict[params['nr_color']]
params = create_params()

# set file_name of hdf5 data file
hdf5file = data_dir + 'Colored_Noise_' + noise_color + '_' + str(params['lowcut']) + '-' + str(params['highcut'])
if simulate and os.path.exists(hdf5file+'.hdf5'):
	print(f'Deleting {hdf5file}.hdf5')
	os.remove(hdf5file+'.hdf5')

# set plot file name
pdffile = plot_dir + 'Colored_Noise_Trial_' + noise_color + '_' + str(params['lowcut']) + '-' + str(params['highcut'])

# set up for loop over run variable
run_range = np.arange(0.8, 1.3, 0.1)
run_variable = 'XL'

for run_value in run_range:
	params[run_variable] = run_value
	run_name = run_variable + '_' + str(run_value).replace('.','')
	hdf5node = run_name

	print(f'running simulation for run variable {run_variable} = {run_value}')
	
	# variable for parallel simulation
	parallel_var = 'noise_level'
	parallel_range=np.arange(0.0025, 0.1, 0.0025) # defines the range of noise level

	# Create an instance of callback class
	nr_simulations = parallel_range.shape[0]
	dc = DataContainer(hdf5file + '.hdf5')

	# Run simulation
	if simulate:
		run_parallel_integration(func, params, parallel_var, parallel_range, hdf5file, hdf5node)
	
	# Create instance of data analyzer
	da = DataAnalyzer(dc)
	noise_freq_range = [params['lowcut'], params['highcut']]
	da.plot_activities(plot_file_name = pdffile + '_' + run_name + '.pdf', nr_variables = params['dimension'], run_name = run_name, sort_variable = parallel_var, noise_freq_range=noise_freq_range)
