limits=( "'noise_highcut':300" "'noise_highcut':400" )
colors=( "'nr_color':1" "'nr_color':2" "'nr_color':3" )


for color in "${colors[@]}"; do
	for limit in "${limits[@]}"; do
		old_limit=$(grep -o \'noise_highcut\':[0-9]* versions/config.py)
		sed -i "s/$old_limit/$limit/g" versions/config.py

		old_color=$(grep -o \'nr_color\':[0-9]* versions/config.py)
		sed -i "s/$old_color/$color/g" versions/config.py

		echo "#!/bin/bash" > jobfile.sh
		echo "#SBATCH --job-name AttractorSimulation" >> jobfile.sh
		echo "#SBATCH --time 2:0:0" >> jobfile.sh
		echo "#SBATCH --cpus 1" >> jobfile.sh
		echo "#SBATCH --mem 2GB" >> jobfile.sh
		echo "#SBATCH --mail-type NONE" >> jobfile.sh
		echo "#SBATCH --output /home/mpib/kamp/LNDG/Noise_Color_Attractor_Model/logs/slurm-%j.out" >> jobfile.sh

		echo "cd $HOME/LNDG/Noise_Color_Attractor_Model" >> jobfile.sh
		echo "module load conda" >> jobfile.sh
		echo "conda activate py3" >> jobfile.sh

		echo "python versions/run_simulation.py" >> jobfile.sh
		sbatch jobfile.sh
		rm jobfile.sh
		sleep 15
	done
done

