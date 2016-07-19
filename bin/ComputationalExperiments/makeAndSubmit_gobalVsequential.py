# Create SLURM scripts and run them on the Rivanna cluster
import os
import time
from datetime import datetime

#-----------------------------------------------------------------
# Check how many jobs are in the queue
#-----------------------------------------------------------------
def checkNumJobsInQueue():
	os.system('squeue -u mb3ad > tmpQueue.txt')
	tmpQueueFile = open('tmpQueue.txt','r')
	count = 0
	for line in tmpQueueFile:
		if line.find('mb3ad') > -1:
			count += 1
	tmpQueueFile.close()
	os.system('rm tmpQueue.txt')
	
	return count
#-----------------------------------------------------------------

numGCs = [2,5,10,15,20,25,30]
numModels = 30

for numgc in numGCs:
	tmpSLURM_name = 'tmpSlurm_' + str(numgc) + '.slurm'
	
	for i in range(numModels):
		tmpSLURM = open(tmpSLURM_name,'w')
		tmpSLURM.write('#!/bin/bash\n')
		tmpSLURM.write('#SBATCH -N 1\n')
		tmpSLURM.write('#SBATCH --ntasks-per-node=16\n')
		tmpSLURM.write('#SBATCH -t 24:00:00\n')
		tmpSLURM.write('#SBATCH -p serial\n')
		tmpSLURM.write('#SBATCH -A CSBLRivanna\n\n')
		
		tmpSLURM.write('module load matlab\n')
		tmpSLURM.write('module load gurobi/6.5.1\n\n')
		
		filePath = '\'/scratch/mb3ad/CE2_globalVsequential_' + str(numgc) + '_gcs' + str(i) + '.tsv\'';
																																						#   N_gcs				rndSeed			resultsFileName				
		tmpSLURM.write('matlab -nosplash -nodesktop -r "addpath(genpath(\'/nv/blue/mb3ad/Desktop/\'));addpath(genpath(\'/share/apps/gurobi/\'));CE2_hpc(' + str(numgc) + ',' + str(i*10) + ',' + filePath + ');quit"\n')
		
		tmpSLURM.close()
		
		# Wait for queue to be less than 50 jobs long
		jobsInQueue = checkNumJobsInQueue()
		print jobsInQueue
		
		while jobsInQueue > 50:
			print 'Waiting...'
			time.sleep(5)
			jobsInQueue = checkNumJobsInQueue()
		
		os.system('sbatch ' + tmpSLURM_name)
		time.sleep(1)
		os.system('rm ' + tmpSLURM_name)

