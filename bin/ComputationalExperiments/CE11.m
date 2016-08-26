% Computational Experiment
% Evaluate iPAU1129 performance predicting gene essentiality (for
% benchmarking purposes)
%
% Written by Matt Biggs, 2016

load PA14_iPAU1129

%## CF medium already defined
modscfm;

% Make gene essentiality predictions
% geneEssentialityiPAU = zeros(size(modPA14_v24.genes));
% for k = 1:length(modPA14_v24.genes)
%     curMod = modPA14_v24;
%     curGene = modPA14_v24.genes{k};
% 
%     delMod = simulateGeneDeletion(curMod,curGene);
%     delGrowth = fba_flex(delMod,modPA14_v24.rxns,modscfm.lb,0);
%     geneEssentialityiPAU(k,1) = delGrowth < 1e-10;
% end
% save('CE11_iPAU_cfGeneEssentiality.mat','geneEssentialityiPAU');

load CE11_iPAU_cfGeneEssentiality.mat

% Load experimental gene essentiality data
fid = fopen('pnas.1419677112.sd03_smaller.txt','r');
genes = textscan(fid, '%s%s%s%s%s%s','Delimiter','\t');
fclose(fid);

geneList = genes{1};
essentialGenes1 = strfind(genes{3},'Essential');
essentialGenesIndicator = ~cellfun(@isempty,essentialGenes1);
essentialGenes = geneList(essentialGenesIndicator);
essentialGenes(1) = [];

essential_genes_in_iPAU = ismember(modPA14_v24.genes,essentialGenes);

% Check accuracy (compared to experimental data)
TP = sum(geneEssentialityiPAU(essential_genes_in_iPAU == 1) == 1);
FP = sum(geneEssentialityiPAU(essential_genes_in_iPAU == 0) == 1);
TN = sum(geneEssentialityiPAU(essential_genes_in_iPAU == 0) == 0);
FN = sum(geneEssentialityiPAU(essential_genes_in_iPAU == 1) == 0);
iPAUaccuracy = (TP + TN)/ (TP + FP + TN + FN)
iPAUprecision = TP / (TP + FP)
iPAUrecall = TP / (TP + FN)


