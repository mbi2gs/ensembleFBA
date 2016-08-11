function [] = CE8_hpc(N_gcs,fractionGenomeAnnotations,rndSeed,resultsFileName)
% Computational Experiment (HPC part)
% Generate networks for CE8
%
% Written by Matt Biggs, 2016

rng(rndSeed,'twister');

% Load universal reaction database and add exchange rxns
load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Get the PA14 data formatted with work with the SEED database
%       PA14Data.biomassFn,growthCarbonSources,growthConditions,nonGrowthCarbonSources,nonGrowthConditions
[PA14Data] = getPA14GrowthConditions(seed_rxns_mat);

% Get the PA14 gene-to-reaction mappings
%       PA14GenomicData.rxn_GPR_mapping
[PA14GenomicData] = getPA14GenomeAnnotations();
rxnList = PA14GenomicData.rxn_GPR_mapping.rxns;
p1 = randperm(length(rxnList), ceil(fractionGenomeAnnotations*length(rxnList)));
rxnList = rxnList(p1);

% Force networks to contain reactions annotated from the genome
Urxns2set = [find(ismember(seed_rxns_mat.rxns,rxnList)); find(ismember(seed_rxns_mat.rxns,'rxn05064'))]; % include spontaneous rxn05064
Uset2 = ones(size(Urxns2set));

% Include exchange reactions for all non-growth conditions, just so that
% it's the network itself--not the lack of exchange reactions--that prevents growth
Xrxns2set = find(sum( abs(seed_rxns_mat.X([PA14Data.growthCarbonSources(:); PA14Data.nonGrowthCarbonSources(:)],:)) ,1) > 0);
Xset2 = ones(size(Xrxns2set));

% Set parameters
biologicalData = struct;
biologicalData.biomassFn = PA14Data.biomassFn;
biologicalData.Urxns2set = Urxns2set;
biologicalData.Uset2 = Uset2;
biologicalData.Xrxns2set = Xrxns2set;
biologicalData.Xset2 = Xset2;

params = struct;
params.sequential = 1;
params.stochast = 1;
params.numModels2gen = 1;
params.verbose = 1;

%------------------------------------------------------------------------
% Gap fill sequentially
%------------------------------------------------------------------------
N_iter = 1;
fprintf(['Number of growth conditions: ' num2str(N_gcs) '\n']);

% Randomly select growth conditions
rp = randperm(size(PA14Data.growthConditions,2),N_gcs);
biologicalData.growthConditions = PA14Data.growthConditions(:,rp);

rp2 = randperm(size(PA14Data.nonGrowthConditions,2),10);
biologicalData.nonGrowthConditions = PA14Data.nonGrowthConditions(:,rp2);

% Sequential gap fill
tic
[modelList_seq] = build_network(seed_rxns_mat,biologicalData,params);
stseq1 = toc;

[gc_bm_vals,ngc_bm_vals] = testModelsInGrowthConditions_flex(modelList_seq,seed_rxns_mat.Ex_names,PA14Data.growthConditions,PA14Data.nonGrowthConditions,0);

m = modelList_seq{1};
m.solveTime = stseq1;
m.gc_bm_vals = gc_bm_vals;
m.ngc_bm_vals = ngc_bm_vals;

save(resultsFileName,'m');

end
