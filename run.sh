# shell for the job:
#PBS -S /bin/bash
#PBS -lnodes=1 -lwalltime=1:59:00 
# job requires at most 30 hours, 0 minutes
#     and 0 seconds wallclock time
# cd to the directory where the program is to be called:
source $HOME/.bash_profile
cd $HOME/projects/rivalry_modeling
# call the programs
echo "Job $PBS_JOBID started at `date`" | mail $USER -s "Job $PBS_JOBID"

ipython new_minimal.py

wait          # wait until programs are finished

echo "Job $PBS_JOBID finished at `date`" | mail $USER -s "Job $PBS_JOBID"
