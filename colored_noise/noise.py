def load_noise(params):
    import numpy as np
    noise_dict = {1:'white', 2:'pink', 3:'blue'}
    noise_color = noise_dict[params['nr_color']]
    noise_lowcut = params['lowcut']
    noise_highcut = params['highcut']
    noise_file = 'colored_noise/%s_noise_%d-%d.csv' % (noise_color, noise_lowcut, noise_highcut)
    noise = np.genfromtxt(noise_file, delimiter=',')
    return noise 

def generate_noise():
    #TODO translate matlab generator into python
    pass
