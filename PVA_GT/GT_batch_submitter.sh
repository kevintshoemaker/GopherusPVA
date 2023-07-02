# loop over all run...sh files in the directory, print the filename and submit the jobs to SLURM
#
mkdir submission_scripts
cd submission_scripts
rm *.sh* 

for i in {1..240..1}; do
    cat > run_PVA_GT_${i}.sh << EOT
#!/bin/bash
#SBATCH --job-name=PVA
#SBATCH --mail-user=kjloope@vt.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=/projects/birdnet/PVA/Rcode/logs/PVAmulti_%j.log
#SBATCH --ntasks=128
#SBATCH --mem=216gb
#SBATCH --time=03:00:00
#SBATCH --account="birdnet"
#SBATCH --partition="normal_q"

if [ -z ${HOME+x} ]; then
  export HOME=$(echo ~)
  source /etc/profile
  source /etc/bashrc
  source $HOME/.bashrc
fi

arrayID=$i
cd /projects/birdnet/PVA/Rcode/
module load R/4.1.0-foss-2021a
Rscript --vanilla GT_RunSimulations.R \$arrayID \$run_name


EOT
done
chmod 777 *.sh
##make a list of the files in order starting with 1
FILES=($(ls -v run*.sh))
##run them one at a time
for FILE in ${FILES[@]}; do
    echo ${FILE}
    echo $1
    sbatch --export=run_name=$1 ${FILE}
    sleep 1 # pause to be kind to the scheduler
done


