#READ ME FIRST
## Submit this job using "qsub submit_job.pbs"
## Queue it will run in
#PBS -q little
## Select 1 node with 1 processor
#PBS -l select=1:ncpus=1
## Pack all of them in 1 node
#PBS -l place=pack
## Name of the Job
#PBS -N test
## Join output and error in a single file
#PBS -j oe
## Export the environment vaiables from your shell
#PBS -V
cd $PBS_O_WORKDIR
## Comment whichever is not required
## Job with 40 MPI processes
#aprun -n 40 your_mpi_executable
## Job with 1 MPI Process and 1 thread
./a.out
