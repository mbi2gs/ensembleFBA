# Generate a list of genes for each species which have close matches to drug targets
#
# Written by Matt Biggs, 2016

import re

#--------------------------------------
# Identify target/drug pairs
#--------------------------------------
targetDictionary = {}
dfile = open('db_small_molecule_targets_protein.fasta','r')
for line in dfile:
	if line.find('>') > -1:
		# Get target ID
		lps = line.split()
		targetID = lps[0].lstrip('>')
		
		# Get drug list
		lps2 = line.split('(DB')
		uflist = 'DB' + lps2[1]
		dlist = uflist.replace(')','').replace('; ',',').rstrip()
		
		# If target is duplicate, combined drug lists
		if targetID in targetDictionary.keys():
			tmpdlist = targetDictionary[targetID].split(',')
			tmpdlist2 = dlist.split(',')
			combinedList = tmpdlist + tmpdlist2
			combinedList = list(set(combinedList))
			targetDictionary[targetID] = ','.join(combinedList)
		else:
			targetDictionary[targetID] = dlist
dfile.close()
#--------------------------------------

print len(targetDictionary)

resultsFiles = ['s.vestibularis.results','s.pneumoniae.results','s.oralis.results','s.mitis.results','s.gallolyticus.results','s.equinus.results']

for rf in resultsFiles:
	# Find all targets (some redundancy)
	file = open(rf,'r')
	targetList = {}
	for line in file:
		lps = line.split()
		m = re.search('peg.\d+',lps[0])
		drugsThatInteractWithTarget = targetDictionary[lps[1]]
		
		# Store target/drug matches
		if m.group(0) in targetList.keys():
			tmpdlist = targetList[m.group(0)].split(',')
			tmpdlist2 = drugsThatInteractWithTarget.split(',')
			combinedList = tmpdlist + tmpdlist2
			combinedList = list(set(combinedList))
			targetList[m.group(0)] = ','.join(combinedList)
		else:
			targetList[m.group(0)] = drugsThatInteractWithTarget
	file.close()
	
	# Record target/drug matches	
	outfile = open(rf + '.targets','w')
	for k in targetList.keys():
		outfile.write( k + '\t' + targetList[k] + '\n')
	outfile.close()
	
	
	
	
	
	
	