function [rxnEssentiality,nonExRxns] = calcStrepEssentialRxns_hpc(modelID, outfileName)
%----------------------------------------------------------------------
% calcStrepEssentialRxns_hpc - Calculates the essential reactions in
% modelID by running FBA seqentially on each one.
%
% Inputs:
%   modelID = Name of mat file containing the current model
%   outfileName = Name of output file where essential rxns will be written
% Outputs:
%   Writes the essential reactions to a text file
%
% Written by Matt Biggs, 2016
%----------------------------------------------------------------------

% Load universal reaction database and add exchange rxns
load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Get the Streptococcus data formatted to work with the SEED database
[StrepData] = getStrepGrowthConditions(seed_rxns_mat);

% Load the current model
eval(['load ' modelID]);

% Create rich media representation
richMedia = sum(StrepData.growthConditions,2);
richMedia(richMedia < -1000) = -1000;

% Simulate all reaction KOs and record essentiality
nonExRxns = m.rxns(~ismember(m.rxns,seed_rxns_mat.Ex_names));
nonExRxns = nonExRxns(~ismember(nonExRxns,'Biomass'));
rxnEssentiality = zeros(size(nonExRxns));
for k = 1:length(nonExRxns)
    curMod = m;
    curRxn = nonExRxns{k};
    curRxnIndex = find(ismember(m.rxns,curRxn));

    curMod.S(:,curRxnIndex) = [];
    curMod.rxns(curRxnIndex) = [];
    curMod.rxnNames(curRxnIndex) = [];
    curMod.rev(curRxnIndex) = [];
    curMod.c(curRxnIndex) = [];
    curMod.ub(curRxnIndex) = [];
    curMod.lb(curRxnIndex) = [];
    curMod.grRules(curRxnIndex) = [];
    
    delGrowth = fba_flex(curMod,seed_rxns_mat.Ex_names,richMedia,0);
    rxnEssentiality(k) = delGrowth < 1e-10;
end

% Write results to file
fid = fopen(outfileName,'w');
for i = 1:length(nonExRxns)
    fprintf(fid,[nonExRxns{i} '\t' num2str(rxnEssentiality(i)) '\n']);
end
fclose(fid);

end