ranges=( "noise_highcut=300" "noise_highcut=400" "noise_highcut=495" )
color=( "colored_noise_gaba_trial.py" )


for color in "${colors[@]}"; do
	for range in "${ranges[@]}"; do
		old_range=$(grep "noise_highcut=" model_config.py)
		sed -i "s/$old_range/$range/g" model_config.py

		old_color=$(grep "noise_color=" model_config.py)
		sed -i "s/$old_color/$color/g" model_config.py

		echo "#!/bin/bash" > jobfile.sh
		echo "#SBATCH --job-name AttractorSimulation" >> jobfile.sh
		echo "#SBATCH --time 2:0:0" >> jobfile.sh
		echo "#SBATCH --cpus 1" >> jobfile.sh
		echo "#SBATCH --mem 2GB" >> jobfile.sh
		echo "#SBATCH --mail-type NONE" >> jobfile.sh
		echo "#SBATCH --output /home/mpib/kamp/LNDG/Noise_Color_Attractor_Model/logs/slurm-%j.out" >> jobfile.sh

		echo "cd $HOME/LNDG/Noise_Color_Attractor_Model" >> jobfile.sh
		echo "module load conda" >> jobfile.sh
		echo "conda activate py2" >> jobfile.sh

		echo "ipython versions/colored_noise_gaba_trial.py" >> jobfile.sh
		sbatch jobfile.sh
		rm jobfile.sh
		sleep 20
	done
done

