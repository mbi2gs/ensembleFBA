% Computational Experiment
% Determine the effect of ensemble size on performance
%
% Written by Matt Biggs, 2016

% Load ensemble
if ~exist('CE7_ensemble_25gcs','var')
    load CE7_ensemble_25gcs
    CE7_ensemble_25gcs(cellfun(@isempty,CE7_ensemble_25gcs)) = [];
end

% Get the PA14 gene-to-reaction mappings
%       PA14GenomicData.rxn_GPR_mapping
[PA14GenomicData] = getPA14GenomeAnnotations();
rxnList = PA14GenomicData.rxn_GPR_mapping.rxns;

% Add GPRs to models in the ensemble
[CE7_ensemble] = addGPRs(CE7_ensemble_25gcs,PA14GenomicData.rxn_GPR_mapping);

% Get complete list of genes from models
N = length(CE7_ensemble);
allGenes = cell(0,1);
for i = 1:N
    allGenes = [allGenes; CE7_ensemble{i}.genes];
end
allGenes = unique(allGenes);

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

%-------------------------------------------------------
% Evaluate ensemble precision/recall
%-------------------------------------------------------
% 1: "Any" threshold
% 2: "Majority" threshold
% 3: "Consensus" threshold

ensembleSizes = 3:2:51;
samples = 10000;
for t = 1:3 % Iterate through each threshold
    ensemble_Accuracy = zeros(length(ensembleSizes),samples);
    ensemble_Precision = zeros(length(ensembleSizes),samples);
    ensemble_Recall = zeros(length(ensembleSizes),samples);
    
    for i = 1:length(ensembleSizes)
        % Choose the threshold
        threshold = 0;
        if t == 1
            threshold = 1;
        elseif t == 2
            threshold = ensembleSizes(i)/2;
        elseif t == 3
            threshold = ensembleSizes(i);
        end
        
        [t threshold]
        
        % Calculate accuracy, precision and recall
        for j = 1:samples
            subset = datasample(1:51,ensembleSizes(i));
            ensembleSum = sum(geneEssentialityByNet(:,subset),2) >= threshold;
            TP = sum( ensembleSum(essential_pegs_inModels == 1) ==  1);
            TN = sum( ensembleSum(essential_pegs_inModels == 0) ==  0);
            FP = sum( ensembleSum(essential_pegs_inModels == 0) ==  1);
            FN = sum( ensembleSum(essential_pegs_inModels == 1) ==  0);
            ensemble_Accuracy(i,j) = (TP + TN) / (TP + TN + FP + FN); % accuracy
            ensemble_Precision(i,j) = TP / (TP + FP); % precision
            ensemble_Recall(i,j) = TP / (TP + FN); % recall
        end
    end

    ensemble_Accuracy = [ensembleSizes(:) ensemble_Accuracy];
    ensemble_Precision = [ensembleSizes(:) ensemble_Precision];
    ensemble_Recall = [ensembleSizes(:) ensemble_Recall];

    % Write to file
    if t == 1
        dlmwrite('CE9_ensembleAccuracyBySize_any.tsv',ensemble_Accuracy,'\t');
        dlmwrite('CE9_ensemblePrecisionBySize_any.tsv',ensemble_Precision,'\t');
        dlmwrite('CE9_ensembleRecallBySize_any.tsv',ensemble_Recall,'\t');
    elseif t == 2
        dlmwrite('CE9_ensembleAccuracyBySize_majority.tsv',ensemble_Accuracy,'\t');
        dlmwrite('CE9_ensemblePrecisionBySize_majority.tsv',ensemble_Precision,'\t');
        dlmwrite('CE9_ensembleRecallBySize_majority.tsv',ensemble_Recall,'\t');
    elseif t == 3
        dlmwrite('CE9_ensembleAccuracyBySize_consensus.tsv',ensemble_Accuracy,'\t');
        dlmwrite('CE9_ensemblePrecisionBySize_consensus.tsv',ensemble_Precision,'\t');
        dlmwrite('CE9_ensembleRecallBySize_consensus.tsv',ensemble_Recall,'\t');
    end
end


