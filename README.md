# Noise Color Attractor Model

Python project to investigate the effect of differently colored noise on a dynamical systems which simulate bistable perception.

New models can be implemented in the versions folder following the syntax of the existing 
model. We simulate the effect of noise on the model by adding noise of different colors to each timestep.

The config.py file configures the model, the parameter of the model. The parameter grid runs serially over all possible parameter settings.  
The parameter settings of the parallel variable are run in parallel. 

The job_wrapper.sh file is a shell file that runs the python scripts through the Slurm batch system, it parallelizes processes over different 
noise color and noise frequency cut of settings.

## Requirements
PyTables, PyGSL, Joblib, scikit-learn, numpy, scipy, matplotlib, seaborn
The repository contains a yml file to set up a conda environment.
