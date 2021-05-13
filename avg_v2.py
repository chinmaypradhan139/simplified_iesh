# Time-Averaging over trajectories run on different cores
import numpy as np
import matplotlib.pyplot as plt 

def read_file(fpath):

	# Read the file
	input_file = open(fpath, 'r') 			# open the file
	lines = input_file.readlines()          # read the content
	input_file.close()                      # close the file
	# print(lines)

	# Store the input in list
	line_list= []
	for line in lines:
	    v_line = line.strip()
	    if (len(v_line) > 0):               # Check and remove blank lines
	        line_list.append(v_line.split())  
	
	print(len(line_list))		# length = 5000
	line_zip = [np.array(list(a), dtype=np.float) for a in  zip(*line_list)] 
	return line_zip

def col_stack(a,b):
	return np.column_stack((a,b))

def time_avg(p,w,N):
	t_avg = np.sum(p,axis=1)*(w/N)#(ntraj/ttraj)
	return t_avg

def write_file(t_p,t_avg):
	# Writing to a file
	f = open("output.txt", "a")
	for i in range(0,tsteps):
		f.write(str(t_p[0][i])+ "   "+str(t_avg[i])+"\n")
	f.close()

def plot_pp(t_p,t_avg):

	plt.plot(t_p[0], t_avg, label = "LUMO population", linestyle = "-", marker = "o")

	plt.xlabel('time')
	plt.ylabel('population')
	plt.legend()
	plt.show()
	plt.savefig('LUMO_pp.png')


#===========================================================

# Read ./i/input.txt
num_core =20
ntraj    = 4
ttraj = 80
tsteps = 90000
#ttraj, ntraj ; dt, ttime;
pp = np.zeros(tsteps)

for seed in range(1,num_core+1):
	fpath = './' +str(seed) +'/fort.21'
	t_p = (read_file(fpath))
	pp = col_stack( pp, t_p[1] )

avg_pp = time_avg(pp, ntraj, ttraj) 
write_file(t_p,avg_pp)	
# plot_pp(t_p,avg_pp)
#============================================================
'''
Test
a = [1,2,3,4,5]
b = [10,20,40,50,100]
c = [8,9,3,5,4]

print(a)
print(b)
print(c)

print(col_stack(a,b))
pp = col_stack(col_stack(a,b),c)
print(pp)

print(time_avg(pp,10,30))
'''
