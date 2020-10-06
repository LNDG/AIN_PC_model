def testing(**kwargs):
	d = {'a':1}
	for key, value in kwargs.items():
		d.update({key:value})
	return d

def params_generator(params, variable, variable_range):
	"""""
	Generator to update model parameters in each iteration
	"""""
	assert (variable in params.keys()), 'Variable not found.'	
	for v in variable_range:
		yield params.update({variable: v})
		 
	
params = {'a':1,'b':2}
import numpy as np
for _ in params_generator(params, 'a', np.arange(10)):
	print(params)
