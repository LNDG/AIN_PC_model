
class model():
	def __init__(self):
		# set default parameters for all models
		self.params = {'nr_timepoints':3000, 'XL': 1.0, 'XR': 1.0, 'gamma': 3.25, 'exponent': 1.0, 'alpha': 1.75, 'tau': 100.0, 'NRa': 2.0, 'NRs': 1.0}
		self.noise_params = {'lowcut': 5, 'highcut': 300, 'nr_color': 1, 'noise_level':0.001}
	
class gaba(model):
	def __init__(self):
		# load and update general parameter settings
		super().__init__()
		self.params.update({'dimension':4}) 
		self.noise_params.update({'nr_noise_tc':2})
		
		from colored_noise.noise import load_noise
		self.noise = load_noise(self.noise_params)

		# initial values
		import pygsl._numobj
		self.init_values = pygsl._numobj.array((0.5, 0.5, 0., 0.))
	
	def func(t, y):
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
		
		dydt[H1] = self.params['XL'] - (1. + y[A1]) * y[H1] - self.params['gamma'] * S(y[H2], self.params['NRa'], self.params['NRs']) 
		dydt[H2] = self.params['XR'] - (1. + y[A2]) * y[H2] - self.params['gamma'] * S(y[H1], self.params['NRa'], self.params['NRs'])
		dydt[A1] = ( -pow(y[A1],self.params['exponent']) + ( self.params['alpha'] * S(y[H1], self.params['NRa'], self.params['NRs']) ) ) / self.params['tau']
		dydt[A2] = ( -pow(y[A2],self.params['exponent']) + ( self.params['alpha'] * S(y[H2], self.params['NRa'], self.params['NRs']) ) ) / self.params['tau']

		return dydt
   
	def create_noise_step(self, iters):
		import numpy as np
		noise_step = np.array([self.noise[0,iters] * self.noise_params['noise_level'], self.noise[1,iters] * self.noise_params['noise_level'], 0.0, 0.0])
