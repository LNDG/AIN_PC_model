import numpy as np
simulate=True
nr_timepoints = 60000
noise_level_range = np.linspace(0.0025, 0.1, 10) # defines the range of noise level
noise_lowcut=10
noise_highcut = 300
mu = {'XL': 1.0, 'XR': 1.0, 'beta': 0.24, 'gamma': 3.75, 'exponent': 1.0, 'alpha': 1.75, 'tau': 100.0, 'NRa': 2.0, 'NRs': 1.0, 'noise_level': 0.0025, 'var_inh': 120.0, 'tau_inh': 50, 'var_inh_infl': 0.0, 'NRa_var_inh': 2.0, 'NRs_var_inh': 1.0, 'var_inh_noise_level': 0.000, 'var_inh_noise_infl': 0.00001, 'noise_color':1}
mu.update({'noise_lowcut': noise_lowcut, 'noise_highcut': noise_highcut})
