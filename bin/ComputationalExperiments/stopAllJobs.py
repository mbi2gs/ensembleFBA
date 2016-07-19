# Create SLURM scripts and run them on the Rivanna cluster
import os
import time
from datetime import datetime

#-----------------------------------------------------------------
# Stop all jobs in the queue
#-----------------------------------------------------------------
os.system('squeue -u mb3ad > tmpQueue.txt')
tmpQueueFile = open('tmpQueue.txt','r')
count = 0
for line in tmpQueueFile:
	if line.find('mb3ad') > -1:
		count += 1
		lps = line.split()
		os.system('scancel ' + lps[0])
tmpQueueFile.close()
os.system('rm tmpQueue.txt')

return count
