modelfiles=( "colored_noise.py" "colored_noise_gaba.py" "colored_noise_feedback.py" )

for model in "${modelfiles[@]}"; do 
	echo "#!/bin/bash" > jobfile.sh
	echo "#SBATCH --job-name AttractorSimulation" >> jobfile.sh
	echo "#SBATCH --time 2:0:0" >> jobfile.sh
	echo "#SBATCH --cpus 5" >> jobfile.sh
	echo "#SBATCH --mem 2GB" >> jobfile.sh
	echo "#SBATCH --mail-type NONE" >> jobfile.sh
	echo "#SBATCH --output /home/mpib/kamp/LNDG/AttractorModel/logs/slurm-%j.out" >> jobfile.sh

	echo "cd $HOME/LNDG/AttractorModel/ppattractor" >> jobfile.sh
	echo "module load conda" >> jobfile.sh
	echo "conda activate py2attractor" >> jobfile.sh

	echo "ipython versions/$model" >> jobfile.sh
	sbatch jobfile.sh
	rm jobfile.sh
done

