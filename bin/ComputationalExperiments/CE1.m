% Computational Experiment
% Show that order matters when gap filling sequentially
%
% Written by Matt Biggs, 2016

% Load universal reaction database and add exchange rxns
load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Get the PA14 data formatted with work with the SEED database
[PA14Data] = getPA14GrowthConditions(seed_rxns_mat);

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
% Does order matter? 
% Gap fill sequentially, in different orders
%------------------------------------------------------------------------
N_iter = 30;
N_gcs_list = [2,5,10,15,20,25,30];
for k = 1:length(N_gcs_list);
    N_gcs = N_gcs_list(k);
    fprintf(['Number of growth conditions: ' num2str(N_gcs) '\n']);
    sims = zeros(N_iter,6); % Jaccard dist, #average unique rxns, # rxns in sequential, # rxns in global, solve time seq (sec), solve time global (sec)
    for i = 1:N_iter
        fprintf(['\titeration ' num2str(i) '\t']);
        % Randomly select growth conditions and 2 permutations
        rp = randperm(size(PA14Data.growthConditions,2),N_gcs);
        randomGrowthConditions = PA14Data.growthConditions(:,rp);
        [p1,p2] = uniquePerms(N_gcs);

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
        
         % Jaccard similarity
        jaccard_sim = jaccardSim(modelList_p1{1}.rxns,modelList_p2{1}.rxns);
        uniqueA = sum(~ismember(modelList_p1{1}.rxns,modelList_p2{1}.rxns));
        uniqueB = sum(~ismember(modelList_p2{1}.rxns,modelList_p1{1}.rxns));
        fprintf(['\tJaccard sim = ' num2str(jaccard_sim) '\n']);

        sims(i,1) = jaccard_sim;
        sims(i,2) = mean([uniqueA,uniqueB]);
        sims(i,3) = length(modelList_p1{1}.rxns);
        sims(i,4) = length(modelList_p2{1}.rxns);
        sims(i,5) = stseq1;
        sims(i,6) = stseq2;
    end
    dlmwrite(['CE1_gapFillSequence_' num2str(N_gcs) '_gcs.tsv'], sims, '\t');
end


