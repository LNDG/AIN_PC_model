class noise(): 
    def __init__(self, params): 
        import numpy as np
        self.params = params
        noise_dict = {1:'white', 2:'pink', 3:'blue'}
        noise_color = noise_dict[params['nr_color']]
        noise_lowcut = params['lowcut']
        noise_highcut = params['highcut']
        self.csvfile = 'colored_noise/%s_noise_%d-%d.csv' % (noise_color, noise_lowcut, noise_highcut)
        self.traces = self.load() 
    
    def load(self):
        noise = np.genfromtxt(self.csvfile, delimiter=',')
        return noise

    def generate(self):
        #TODO translate matlab generator into python
        pass

    def create_step(self, iters):
        assert self.traces.shape[0] > 2, 'Noise must have more than 2 rows.'
        step = np.zeros(self.params['nr_noise_tc'] + 2) # plus 2 because no noise is added to adaption variable
        for tc in range(self.params['nr_noise_tc']):
            step[tc] = self.traces[tc, iters] * self.params['noise_level']
        return step


