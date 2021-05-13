import os
import subprocess
# import stat

num_core = 1
tot_traj = 1					# Each core runs (tot_traj/num_core) trajectories
# dt, ttime

compiled_fname = "a.out"
job_script     = "job.sh"
input_filename = "input.txt" 

#Input text files
for seed in range(1,num_core+1):

	file = open(input_filename, 'w')
	file.write( str(num_core)+" "+str(tot_traj)+" "+ str(seed) +"\n")
	file.close()

	file = open("bashrun.sh", "w")
	file.write("#!/bin/bash\n")
	file.write("mkdir "+ str(seed) +"\n")
		#check ./a.out present
	file.write("cp "+ compiled_fname + " " + str(seed) + "\n")
	file.write("cp "+ job_script	  + " " + str(seed) + "\n")
	file.write("cp "+ "fort.23" + " " + str(seed) + "\n")
	file.write("cp "+ "raw_x.txt"+ " " + str(seed) + "\n")
	file.write("cp "+ "raw_w.txt"+ " " + str(seed) + "\n")
	file.write("mv "+ input_filename  + " " + str(seed) + "\n")

	file.write("cd "+ str(seed) + "\n")
	file.write("qsub job.sh	\n") 		# ("./"+ output_filename+ " \n")
	file.write("cd .. \n") 
	file.close()

	# print(os.getcwd())
	# st = os.stat(os.getcwd()+'/bashrun.sh')
	# os.chmod(os.getcwd()+'/bashrun.sh', st.st_mode | stat.S_IEXEC)
	
	subprocess.call("(chmod +x bashrun.sh)", shell = True)
	subprocess.call("(./bashrun.sh)", shell = True)






