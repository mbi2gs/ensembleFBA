function [] = CE2_hpc(N_gcs,rndSeed,resultsFileName)
% Computational Experiment
% Compare parsimony between sequential and global gap filling approaches
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
biologicalData.nonGrowthConditions = [];
biologicalData.biomassFn = PA14Data.biomassFn;
biologicalData.Urxns2set = Urxns2set;
biologicalData.Uset2 = Uset2;
biologicalData.Xrxns2set = Xrxns2set;
biologicalData.Xset2 = Xset2;

params = struct;
params.stochast = 0;
params.numModels2gen = 1;
params.verbose = 0;

rng(rndSeed,'twister');

jaccardSim = @(a,b) sum(ismember(a,b))/length(unique([a(:);b(:)]))';

%------------------------------------------------------------------------
% Does order matter? 
% Gap fill sequentially, in different orders
%------------------------------------------------------------------------
N_iter = 1;
fprintf(['Number of growth conditions: ' num2str(N_gcs) '\n']);
sims = zeros(N_iter,6); % Jaccard dist, #average unique rxns, # rxns in sequential, # rxns in global, solve time seq (sec), solve time global (sec)

% Randomly select growth conditions and 2 permutations
rp = randperm(size(PA14Data.growthConditions,2),N_gcs);
biologicalData.growthConditions = PA14Data.growthConditions(:,rp);

% Sequential gap fill
params.sequential = 1;

tic
[modelList_seq] = build_network(seed_rxns_mat,biologicalData,params);
stseq1 = toc;

% Global gap fill
params.sequential = 0;

tic
[modelList_glob] = build_network(seed_rxns_mat,biologicalData,params);
stseq2 = toc;

 % Jaccard similarity
jaccard_sim = jaccardSim(modelList_seq{1}.rxns,modelList_glob{1}.rxns);
uniqueSeq = sum(~ismember(modelList_seq{1}.rxns,modelList_glob{1}.rxns));
uniqueGlob = sum(~ismember(modelList_glob{1}.rxns,modelList_seq{1}.rxns));
fprintf(['Jaccard sim = ' num2str(jaccard_sim) '\n']);

sims(1,1) = jaccard_sim;
sims(1,2) = mean([uniqueSeq,uniqueGlob]);
sims(1,3) = length(modelList_seq{1}.rxns);
sims(1,4) = length(modelList_glob{1}.rxns);
sims(1,5) = stseq1;
sims(1,6) = stseq2;

dlmwrite(resultsFileName, sims, '\t');


end
