% Computational Experiment
% Evaluate sources of model diversity and their effects on ensemble
% performance
%
% Written by Matt Biggs, 2016

% Load ensemble
ensembleNames = {};
for ngcs = numGCs
    for i = 1:length(fractionGenomeAnnotations)
        fga = fractionGenomeAnnotations{i};
        ensembleName = ['CE8_' num2str(ngcs) 'gcs_10ngcs_' fga 'fracGenAnn_ensemble'];
        
        if ~exist(ensembleName,'var')
            eval(['load ' ensembleName])
        end
        
        command = [ensembleName '(cellfun(@isempty,' ensembleName ')) = [];'];
        eval(command);
        ensembleNames{end+1,1} = ensembleName;
    end
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

% Calculate gene essentiality for each of the 16 ensembles
geneEssentialityByNet_16 = cell(16,1);
ensemble_Accuracy_Precision_Recall_16 = zeros(16,3);
for i = 1:length(ensembleNames)
    eval(['currEnsemble = ' ensembleNames{i}]);
    % Add GPRs to models in the ensemble
    [currEnsemble] = addGPRs(currEnsemble,PA14GenomicData.rxn_GPR_mapping);

    % Get complete list of genes from models
    N = length(currEnsemble);
    allGenes = cell(0,1);
    for k = 1:N
        allGenes = [allGenes; currEnsemble{k}.genes];
    end
    allGenes = unique(allGenes);
    
    % Simulate gene knockouts, run FBA and check for growth within all models 
    geneEssentialityByNet = zeros(length(allGenes),N);
    
    for j = 1:length(currEnsemble)
        [curMod,~] = addExchangeRxns(currEnsemble{j}, {'cpd00011','cpd00060','cpd00065','cpd00069','cpd00084','cpd00156','cpd00209'});

        for k = 1:length(allGenes)
            curGene = allGenes{k};
            delMod = simulateGeneDeletion(curMod,curGene);
            delGrowth = fba_flex(delMod,seed_rxns_mat.Ex_names,cfSputum,0);
            geneEssentialityByNet(k,j) = delGrowth < 1e-10;
        end
    end
    geneEssentialityByNet_16{i,1} = geneEssentialityByNet;

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

    % Evaluate ensemble accuracy/precision/recall
    threshold = N/2;
    ensembleSum = sum(geneEssentialityByNet,2) >= threshold;
    TP = sum( ensembleSum(essential_pegs_inModels == 1) ==  1);
    TN = sum( ensembleSum(essential_pegs_inModels == 0) ==  0);
    FP = sum( ensembleSum(essential_pegs_inModels == 0) ==  1);
    FN = sum( ensembleSum(essential_pegs_inModels == 1) ==  0);
    ensemble_Accuracy_Precision_Recall_16(i,1) = (TP + TN) / (TP + TN + FP + FN); % accuracy
    ensemble_Accuracy_Precision_Recall_16(i,2) = TP / (TP + FP); % precision
    ensemble_Accuracy_Precision_Recall_16(i,3) = TP / (TP + FN); % recall
end

save('CE8_geneEssentiality_AccPrecRec.mat','geneEssentialityByNet_16','ensemble_Accuracy_Precision_Recall_16');

% Write to file
% NOTE: Each row corresponds to the ensembles listed in: "ensembleNames"
dlmwrite('CE8_ensemblePrecisionsNRecalls.tsv',ensemble_Accuracy_Precision_Recall_16,'\t');


