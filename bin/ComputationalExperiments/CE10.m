% Computational Experiment
% Generate three example GENREs based only on the order of gap filling
%
% Written by Matt Biggs, 2016

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
params.stochast = 0;
params.numModels2gen = 1;
params.verbose = 0;

jaccardSim = @(a,b) sum(ismember(a,b))/length(unique([a(:);b(:)]))';

%------------------------------------------------------------------------
% Gap fill sequentially, in different orders
%------------------------------------------------------------------------
N_gcs = 10;

% Randomly select growth conditions and 2 permutations
rng(2041);
rp = randperm(size(PA14Data.growthConditions,2),N_gcs);
randomGrowthConditions = PA14Data.growthConditions(:,rp);
[p1,p2] = uniquePerms(N_gcs);
rng(2208);
[p3,~] = uniquePerms(N_gcs);

% Sequential gap fill using order from the first random permutation
biologicalData.growthConditions = randomGrowthConditions(:,p1);
biologicalData.nonGrowthConditions = [];

tic
[modelList_p1] = build_network(seed_rxns_mat,biologicalData,params);
stseq1 = toc;

% Sequential gap fill using order from the second random permutation
biologicalData.growthConditions = randomGrowthConditions(:,p2);
biologicalData.nonGrowthConditions = [];

tic
[modelList_p2] = build_network(seed_rxns_mat,biologicalData,params);
stseq2 = toc;

% Sequential gap fill using order from the third random permutation
biologicalData.growthConditions = randomGrowthConditions(:,p3);
biologicalData.nonGrowthConditions = [];

tic
[modelList_p3] = build_network(seed_rxns_mat,biologicalData,params);
stseq3 = toc;

 % Jaccard similarity
jaccard_sim = jaccardSim(modelList_p1{1}.rxns,modelList_p2{1}.rxns);
uniqueA = sum(~ismember(modelList_p1{1}.rxns,modelList_p2{1}.rxns));
uniqueB = sum(~ismember(modelList_p2{1}.rxns,modelList_p1{1}.rxns));
fprintf(['\tJaccard sim = ' num2str(jaccard_sim) '\n']);

sims = zeros(1,6); % Jaccard similarity, #average unique rxns, # rxns in network1, # rxns in network2, solve time 1 (sec), solve time 2 (sec)
sims(1,1) = jaccard_sim;
sims(1,2) = mean([uniqueA,uniqueB]);
sims(1,3) = length(modelList_p1{1}.rxns);
sims(1,4) = length(modelList_p2{1}.rxns);
sims(1,5) = stseq1;
sims(1,6) = stseq2;

% Extract info for figure 
uniqueRxnsA = modelList_p1{1}.rxnNames(~ismember(modelList_p1{1}.rxns,modelList_p2{1}.rxns));
uniqueRxnsB = modelList_p2{1}.rxnNames(~ismember(modelList_p2{1}.rxns,modelList_p1{1}.rxns));

[gc_bm_vals1,ngc_bm_vals1] = testModelsInGrowthConditions_flex(modelList_p1,seed_rxns_mat.Ex_names,PA14Data.growthConditions,PA14Data.nonGrowthConditions,0);
[gc_bm_vals2,ngc_bm_vals2] = testModelsInGrowthConditions_flex(modelList_p2,seed_rxns_mat.Ex_names,PA14Data.growthConditions,PA14Data.nonGrowthConditions,0);
[gc_bm_vals3,ngc_bm_vals3] = testModelsInGrowthConditions_flex(modelList_p3,seed_rxns_mat.Ex_names,PA14Data.growthConditions,PA14Data.nonGrowthConditions,0);
