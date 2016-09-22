# Parse KEGG.pathways.tsv
#
# Written by Matt Biggs, 2016

# Read in subsystem assignments for KEGG reactions
infile = open('KEGG.pathways.tsv','r')
rxnSubsysDictionary = {}
for line in infile:
	if line.find('Source	Source ID') == -1:
		lps = line.split('\t')
		subsys = lps[2]
		rxns = lps[4].strip().split('|')
		if len(rxns[0]) > 0:
			for r in rxns:
				if r in rxnSubsysDictionary.keys():
					rxnSubsysDictionary[r] += ('|' + subsys)
				else:
					rxnSubsysDictionary[r] = subsys
infile.close()

# Read in KEGG-to-SEED reaction mappings
infile = open('KEGG.aliases.txt','r')
SEEDrxnSubsysDictionary = {}
for line in infile:
	if line.find('KEGG	default	plantdefault') == -1:
		lps = line.split('\t')
		KEGGid = lps[0]
		SEEDids = lps[1].split('|')
		if KEGGid in rxnSubsysDictionary.keys():
			for sid in SEEDids:
				SEEDrxnSubsysDictionary[sid] = rxnSubsysDictionary[KEGGid]
infile.close()

# Write mapping to file
outfile = open('SEED.rxn.pathways.tsv','w')
for r in SEEDrxnSubsysDictionary.keys():
	outfile.write(r + '\t' + SEEDrxnSubsysDictionary[r] + '\n')
outfile.close()


