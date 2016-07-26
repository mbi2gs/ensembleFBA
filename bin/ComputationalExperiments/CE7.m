% Computational Experiment
% Predict gene essentiality
%
% Written by Matt Biggs, 2016

% Load ensemble
if ~exist('CE7_ensemble_25gcs','var')
    load CE7_ensemble_25gcs
    CE7_ensemble_25gcs(cellfun(@isempty,CE7_ensemble_25gcs)) = [];
end

% Get universal rxn database
load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Get the PA14 gene-to-reaction mappings
%       PA14GenomicData.rxn_GPR_mapping
[PA14GenomicData] = getPA14GenomeAnnotations();
rxnList = PA14GenomicData.rxn_GPR_mapping.rxns;

% Load in CF sputum media condition
[cfSputum] = getCFSputumMedium(seed_rxns_mat);

% Add GPRs to models in the ensemble
[CE7_ensemble] = addGPRs(CE7_ensemble_25gcs,PA14GenomicData.rxn_GPR_mapping);

% Get complete list of genes from models
N = length(CE7_ensemble);
allGenes = cell(0,1);
for i = 1:N
    allGenes = [allGenes; CE7_ensemble{i}.genes];
end
allGenes = unique(allGenes);

% Simulate gene knockouts, run FBA and check for growth within all models 
% geneEssentialityByNet = zeros(length(allGenes),N);
% for j = 1:N
%     [curMod,~] = addExchangeRxns(CE7_ensemble{j}, {'cpd00011','cpd00060','cpd00065','cpd00069','cpd00084','cpd00156','cpd00209'});
% 
%     for k = 1:length(allGenes)
%         curGene = allGenes{k};
%         delMod = simulateGeneDeletion(curMod,curGene);
%         delGrowth = fba_flex(delMod,seed_rxns_mat.Ex_names,cfSputum,0);
%         geneEssentialityByNet(k,j) = delGrowth < 1e-10;
%     end
% end
% save('CE7_geneEssentialityByNet.mat','geneEssentialityByNet');

load CE7_geneEssentialityByNet.mat

% Load experimental gene essentiality data
fid = fopen('pnas.1419677112.sd03_smaller.txt','r');
genes_experimental = textscan(fid, '%s%s%s%s%s%s','Delimiter','\t');
fclose(fid);

geneList = genes_experimental{1};
essentialGenes1 = strfind(genes_experimental{3},'Essential');
essentialGenesIndicator = ~cellfun(@isempty,essentialGenes1);
essentialGenes = geneList(essentialGenesIndicator);
essentialGenes(1) = [];

% Load peg-to-PA14 gene ID mappings
fid = fopen('PA14_FeatureTable_smaller.txt','r');
genes_mapping = textscan(fid, '%s%s%s','Delimiter','\t');
fclose(fid);

trimmed_pegs = genes_mapping{1};
trimmed_pegs = cellfun(@(orig,old,new) strrep(orig,old,new), trimmed_pegs, repmat({'fig|208963.12.'},[length(trimmed_pegs),1]), repmat({''},[length(trimmed_pegs),1]),'UniformOutput',false);
trimmed_pegs(1) = [];

pa14_IDs = genes_mapping{2};
pa14_IDs(1) = [];

essential_pegs = trimmed_pegs(ismember(pa14_IDs,essentialGenes));
essential_pegs_inModels = ismember(allGenes,essential_pegs);

% Evaluate the networks individually and as ensembles
network_Accuracy_Precision_Recall = zeros(N,3);
ensemble_Accuracy_Precision_Recall = zeros(3,3);

% Evaluate individual networks
for j = 1:N
    TP = sum( geneEssentialityByNet(essential_pegs_inModels == 1,j) ==  1);
    TN = sum( geneEssentialityByNet(essential_pegs_inModels == 0,j) ==  0);
    FP = sum( geneEssentialityByNet(essential_pegs_inModels == 0,j) ==  1);
    FN = sum( geneEssentialityByNet(essential_pegs_inModels == 1,j) ==  0);
    network_Accuracy_Precision_Recall(j,1) = (TP + TN) / (TP + TN + FP + FN); % accuracy
    network_Accuracy_Precision_Recall(j,2) = TP / (TP + FP); % precision
    network_Accuracy_Precision_Recall(j,3) = TP / (TP + FN); % recall
end

% Evaluate ensemble precision/recall
%%%%%%%%%%%%%%%%%% Stopped working here!
thresholds = [1,N/2,N];
for t = 1:N
    threshold = thresholds(t);
    testGCs = zeros(length(testGCs),N);
    testNGCs = zeros(length(testNGCs),N);
    for j = 1:N
        testGCs(:,j) = curEnsemble{j}.gc_bm_vals(31:end);
        testNGCs(:,j) = curEnsemble{j}.ngc_bm_vals(1:17);
    end
    testGCs = sum(testGCs > 1e-10,2) >= t;
    testNGCs = sum(testNGCs > 1e-10,2) >= t;
    TP = sum( testGCs );
    TN = sum( testNGCs == 0 );
    FP = sum( testNGCs );
    FN = sum( testGCs == 0 );
    ensemblePrecisionsNRecalls(end+1,1) = TP / (TP + FP); % precision
    ensemblePrecisionsNRecalls(end,2) = TP / (TP + FN); % recall
    ensemblePrecisionsNRecalls(end,3) = t;
    ensemblePrecisionsNRecalls(end,4) = gcs;
end


% Write to file
% dlmwrite('CE7_networkAccuracies.tsv',networkAccuracy,'\t');
% dlmwrite('CE7_ensembleAccuracy.tsv',ensembleAccuracy,'\t');
% dlmwrite('CE7_ensemblePrecisionsNRecalls.tsv',ensemblePrecisionsNRecalls,'\t');
% 
% 
