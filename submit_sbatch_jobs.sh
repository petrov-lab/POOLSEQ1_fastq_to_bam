#!/bin/bash
cmd=$1
samp=$2
cmdName=$3
time=$4
mem=$5
tasks=${6:-1}
threads=${7:-1}
depend=${8:-"none"}
jobscriptDir=${9:-$HOME}

echo \
"#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=$cmdName.$samp
#################  
#a file for job output, you can check job progress
#SBATCH --output=$jobscriptDir/$cmdName.$samp.output
#################
# a file for errors from the job
#SBATCH --error=$jobscriptDir/$cmdName.$samp.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=$time
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB=4GB per CPU
#SBATCH --mem=$mem
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#tasks to run per node; a task is usually mapped to a MPI processes.
# for local parallelism (OpenMP or threads), use --ntasks-per-node=1
#SBATCH --ntasks-per-node=$tasks
#################
#number of threads to use per node - only used if the application itself is multi-threaded
#SBATCH --cpus-per-task=$threads
#########################
# partition to run on: owners may get you kickedoff without notice
#SBATCH --partition=dpetrov,hns,normal,owners
#################
#output directory for temporary + intermediate output/error files; could be big
#SBATCH --chdir=$jobscriptDir
#" > $jobscriptDir/$cmdName.$samp.sbatch

if [ ! $depend = "none" ]; then
	echo "#################
#wait for specified job to finish successfully
#SBATCH --dependency=afterok:$depend" >> $jobscriptDir/$cmdName.$samp.sbatch
fi

echo "
#now run normal batch commands
$cmd " >> $jobscriptDir/$cmdName.$samp.sbatch

echo "sbatch script written: $jobscriptDir/$cmdName.$samp.sbatch"
sid=$(sbatch $jobscriptDir/$cmdName.$samp.sbatch | cut -f4 -d' ')
echo $sid
