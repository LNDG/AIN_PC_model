# AIN_PC_model

Python project that implements a dynamical systems model to simulate bistable perception.

Different versions of the model implement different ideas regarding transition-related activations. The over-arching idea is that these models implement specific transition-related activations that change the workings of the model around the time of the transition. In this manner, the model can account for some recent findings between the strength of transition-related responses and perception.  

Results are saved in an hdf5 file, and simulations across parameter settings are run in parallel for speed, within a single node. The run.sh file is a shell file that runs the python scripts through the TORQUE batch system, for parallelization across nodes on a cluster.

## Requirements
PyTables, PyGSL, Parallel Python (pp), numpy, scipy, matplotlib
