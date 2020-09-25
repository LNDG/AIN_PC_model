ranges=( "noise_lowcut=5" "noise_lowcut=10" "noise_lowcut=50" "noise_lowcut=100" )
modelfiles=( "colored_noise_gaba.py" )


for model in "${modelfiles[@]}"; do
	for range in "${ranges[@]}"; do
		old_range=$(grep "noise_lowcut=" model_config.py)
		sed -i "s/$old_range/$range/g" model_config.py
		echo "#!/bin/bash" > jobfile.sh
		echo "#SBATCH --job-name AttractorSimulation" >> jobfile.sh
		echo "#SBATCH --time 2:0:0" >> jobfile.sh
		echo "#SBATCH --cpus 5" >> jobfile.sh
		echo "#SBATCH --mem 5GB" >> jobfile.sh
		echo "#SBATCH --mail-type NONE" >> jobfile.sh
		echo "#SBATCH --output /home/mpib/kamp/LNDG/Noise_Color_Attractor_Model/logs/slurm-%j.out" >> jobfile.sh

		echo "cd $HOME/LNDG/Noise_Color_Attractor_Model" >> jobfile.sh
		echo "module load conda" >> jobfile.sh
		echo "conda activate py2" >> jobfile.sh

		echo "ipython versions/$model" >> jobfile.sh
		sbatch jobfile.sh
		rm jobfile.sh
		sleep 20
	done
done

