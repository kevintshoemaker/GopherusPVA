# loop over all run...sh files in the directory, print the filename and submit the jobs to SLURM
#
mkdir submission_scripts/ 
cd submission_scripts/
rm *.sh* #delete older files
for i in {1..24..1}; do
#make the appropriate submission files.  edit the SBATCH parameters according to your local cluster.
    cat > run_PVA_DT_${i}.sh << EOT
#!/bin/bash
#SBATCH --job-name=PVA_DT
#SBATCH --mail-user=kjloope@vt.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output=/projects/birdnet/PVA_DT/Rcode_DT/logs/PVA_DT_%j.log
#SBATCH --ntasks=128
#SBATCH --mem=216gb
#SBATCH --time=10:00:00
#SBATCH --account="birdnet"
##SBATCH --partition="normal_q"

if [ -z ${HOME+x} ]; then
  export HOME=$(echo ~)
  source /etc/profile
  source /etc/bashrc
  source $HOME/.bashrc
fi

##mkdir -p /fastscratch/kjloope/\$run_name/
arrayID=$i
cd /projects/birdnet/PVA_DT/Rcode_DT/
module load R/4.1.0-foss-2021a
Rscript --vanilla DT_RunSimulations.R \$arrayID \$run_name

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
    sleep 2 # pause to be kind to the scheduler
done


