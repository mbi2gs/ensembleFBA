% Computational Experiment
% Compare parsimony between sequential and global gap filling approaches
%
% Written by Matt Biggs, 2016

% Load universal reaction database and add exchange rxns
load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Get the PA14 data formatted with work with the SEED database
%   Includes: modPA14_v24, modglc, modlb, modscfm
load PA14_iPAU1129

% Make "universal" reaction database that includes the unique reactions 
% from iPAU1129
srm_mets = cellstr([char(seed_rxns_mat.mets) repmat('[c]',size(seed_rxns_mat.mets))]);
allMets = unique([modPA14_v24.mets; srm_mets]);

seed_plus_iPAU = struct;
seed_plus_iPAU.mets = allMets;

% Match metabolite indices between models and databases
iPAU_2_spi = zeros(size(modPA14_v24.mets));
seed_2_spi = zeros(size(srm_mets));
for i = 1:length(seed_plus_iPAU.mets)
   curMet = seed_plus_iPAU.mets{i};
   iPAU_2_spi(ismember(modPA14_v24.mets,curMet)) = i;
   seed_2_spi(ismember(srm_mets,curMet)) = i;
end

% Update metabolite formulas and names
seed_plus_iPAU.metFormulas = cell(size(allMets));
seed_plus_iPAU.metFormulas(iPAU_2_spi) = modPA14_v24.metFormulas;
seed_plus_iPAU.metFormulas(seed_2_spi) = seed_rxns_mat.metFormulas;
seed_plus_iPAU.metNames = cell(size(allMets));
seed_plus_iPAU.metNames(iPAU_2_spi) = modPA14_v24.metNames;
seed_plus_iPAU.metNames(seed_2_spi) = seed_rxns_mat.metNames;

% Add reactions from SEED
seed_plus_iPAU.rxns = seed_rxns_mat.rxns;
seed_plus_iPAU.S = sparse(zeros(length(seed_plus_iPAU.mets),length(seed_plus_iPAU.rxns)));
for i = 1:length(seed_plus_iPAU.rxns)
   seed_plus_iPAU.S(seed_2_spi,i) = seed_rxns_mat.S(:,i);
end

% Add reactions from iPAU1129
unique_iPAU_rxns = modPA14_v24.rxns(~ismember(modPA14_v24.rxns,seed_rxns_mat.rxns));
iPAU_ExchangeRxns = unique_iPAU_rxns(~cellfun(@isempty,strfind(unique_iPAU_rxns,'EX_')));
unique_iPAU_rxns = unique_iPAU_rxns(cellfun(@isempty,strfind(unique_iPAU_rxns,'EX_')));
unique_iPAU_rxns = unique_iPAU_rxns(~ismember(unique_iPAU_rxns,'PA14_Biomass'));

seed_plus_iPAU.rxns = [seed_plus_iPAU.rxns; unique_iPAU_rxns];
seed_plus_iPAU.S = [seed_plus_iPAU.S sparse(zeros(length(seed_plus_iPAU.mets),length(unique_iPAU_rxns)))];
for i = 1:length(unique_iPAU_rxns)
    curRxn = unique_iPAU_rxns{i};
    seed_plus_iPAU.S(iPAU_2_spi,length(seed_rxns_mat.rxns)+i) = modPA14_v24.S(:,ismember(modPA14_v24.rxns,curRxn));
end

% Update rxnNames and reversibility indicators
iPAU_urevs = zeros(size(unique_iPAU_rxns));
iPAU_urxnNames = cell(size(unique_iPAU_rxns));
for i = 1:length(unique_iPAU_rxns)
    curRxn = unique_iPAU_rxns{i};
    iPAU_urevs(i) = modPA14_v24.rev(ismember(modPA14_v24.rxns,curRxn));
    iPAU_urxnNames{i} = modPA14_v24.rxnNames(ismember(modPA14_v24.rxns,curRxn));
end

seed_plus_iPAU.rev = [seed_rxns_mat.rev; iPAU_urevs];
seed_plus_iPAU.rxnNames = [seed_rxns_mat.rxnNames; iPAU_urxnNames];

% Build exchange reaction matrix
seed_plus_iPAU.Ex_names = [seed_rxns_mat.Ex_names; iPAU_ExchangeRxns];
seed_plus_iPAU.X = sparse(zeros(length(seed_plus_iPAU.mets),length(seed_plus_iPAU.Ex_names)));

for i = 1:length(seed_rxns_mat.Ex_names)
   seed_plus_iPAU.X(seed_2_spi,i) = seed_rxns_mat.X(:,i);
end

for i = 1:length(iPAU_ExchangeRxns)
    curRxn = iPAU_ExchangeRxns{i};
    seed_plus_iPAU.X(iPAU_2_spi,length(seed_rxns_mat.Ex_names)+i) = modPA14_v24.S(:,ismember(modPA14_v24.rxns,curRxn));
end

iPAU_biomassRxn = sparse(zeros(length(seed_plus_iPAU.mets),1));
iPAU_biomassRxn(iPAU_2_spi,1) = modPA14_v24.S(:,ismember(modPA14_v24.rxns,'PA14_Biomass'));

% Convert growth conditions to new metabolite list format
% Cobalt                cpd00149    9802    (Index)
% Copper                cpd00058    9803
% Fe3                   cpd10516    9831
% H+                    cpd00067    9850
% H2O                   cpd00001    9851
% Mg2                   cpd00254    9894
% N2                    cpd00528    9904
% NH4                   cpd00013    9900
% O2                    cpd00007    9905
% Pi                    cpd00009    9908
% SO4                   cpd00048    9917
% Cd2+                  cpd01012    9796
% Fe2                   cpd00021    9830
% K+                    cpd00205    9858
% Mn2                   cpd00030    9895
% Na                    cpd00971    9898
% Zn2                   cpd00034    9932
% Growth Carbon sources:
% 4-Hydroxybenzoate     cpd00136    1346
% Acetate               cpd00029    9776
% Citrate               cpd00137    9799
% Glycerol-3-phosphate  cpd00080    9839
% D-Alanine             cpd00117    9809
% D-Fructose            cpd00082    9815
% Adenosine             cpd00182    9778
% Fumarate              cpd00106    9833
% L-Arginine            cpd00051    9863
% L-Asparagine          cpd00132    9864 %
% L-Aspartate           cpd00041    9865
% L-Glutamate           cpd00023    9868
% L-Glutamine           cpd00053    9869
% L-Histidine           cpd00119    9870
% Glycerol              cpd00100    9838 %
% Glycine               cpd00033    9840
% Ornithine             cpd00064    9906
% L-Phenylalanine       cpd00066    9880
% L-Proline             cpd00129    9881
% L-Serine              cpd00054    9882
% Malonate              cpd00308    9888
% L-Alanine             cpd00035    9860
% Glucose               cpd00027    9820
% Putrescine            cpd00118    9912
% Pyruvate              cpd00020    9913
% Succinate             cpd00036    9916
% 2-Oxoglutarate        cpd00024    9773
% GABA                  cpd00281    9938
% L-Isoleucine          cpd00322    9872
% L-Malate              cpd00130    9876
% L-Lactate             cpd00159    9873
% L-Leucine             cpd00107    9874
% D-Mannitol            cpd00314    9935
% Aminoethanol          cpd00162    9786
% N-Acetyl-L-glutamate  cpd00477    3940
% gly-glu-L             cpd11592    9844
% gly-pro-L             cpd11588    9848
% Hydroxy-L-proline     cpd00851    9923
% Butyrate              cpd00211    9793
% L-alanylglycine       cpd11585    9861
% L-Lysine              cpd00039    9875
% Carnitine             cpd00266    6547
% GLCN                  cpd00222    7244
% Hydroxyphenylacetate  cpd00489    9774
% Propionate            cpd00141    9911
% Uridine               cpd00249    9926
% Hydroxybutanoate      cpd00797    9770

minimalMediaBase = zeros(length(seed_plus_iPAU.Ex_names),1);
minimalMediaBase([9802,9803,9831,9850,9851,9894,9904,9900,9905, ...
				  9908,9917,9796,9830,9858,9895,9898,9932],1) = -1000;

growthCarbonSources = [1346,9776,9799,9839,9809,9815,9778,9833,9863, ...
					   9864,9865,9868,9869,9870,9838,9840,9906,9880, ...
					   9881,9882,9888,9860,9820,9912,9913,9916,9773, ...
					   9938,9872,9876,9873,9874,9935,9786,3940,9844, ...
					   9848,9923,9793,9861,9875,6547,7244,9774,9911, ...
					   9926,9770];
                   
n = length(growthCarbonSources);
growthConditions = repmat(minimalMediaBase,[1,n]);
for i = 1:n
    growthConditions(growthCarbonSources(i),i) = -10;
end

% [newBiomassFn] = removeBlockedBiomassComponents(seed_plus_iPAU,growthConditions(:,1),iPAU_biomassRxn);
newBiomassFn = iPAU_biomassRxn;
newBiomassFn([3 250 6246]) = 0;
% Removed: Ubiquinol-9                           cJB00125[c]
%          L-Valine                              cpd00156[c]
%          Peptidoglycan polymer (n subunits)    cpd15665[c]

% Include exchange reactions for all growth conditions
Xrxns2set = growthCarbonSources;
Xset2 = ones(size(Xrxns2set));

iPAU_nonExchangeRxns = modPA14_v24.rxns(cellfun(@isempty,strfind(modPA14_v24.rxns,'EX_')));

% Set parameters
biologicalData = struct;
biologicalData.biomassFn = newBiomassFn;
biologicalData.nonGrowthConditions = [];
biologicalData.Xrxns2set = Xrxns2set;
biologicalData.Xset2 = Xset2;

params = struct;
params.stochast = 0;
params.numModels2gen = 1;
params.verbose = 0;

jaccardSim = @(a,b) sum(ismember(a,b))/length(unique([a(:);b(:)]))';

%------------------------------------------------------------------------
% Does global gap filling produce more biologically relevant networks?
% Gap fill sequentially, in different orders
%------------------------------------------------------------------------
N_iter = 30;
N_gcs = 10;
N_rxns = ceil(length(iPAU_nonExchangeRxns)*0.8);
sims = zeros(N_iter,6);
for i = 1:N_iter
    fprintf(['\titeration ' num2str(i) '\t']);
    % Randomly select growth conditions
    rp = randperm(size(growthConditions,2),N_gcs);
    biologicalData.growthConditions = growthConditions(:,rp);
    
    numericalIssues = 1;
    while numericalIssues == 1
        % Force networks to contain a subset of PA14 reactions
        rxnSubset = randperm(length(iPAU_nonExchangeRxns),N_rxns);
        tmpRxnList = iPAU_nonExchangeRxns(rxnSubset);
        biologicalData.Urxns2set = find(ismember(seed_plus_iPAU.rxns,tmpRxnList));
        biologicalData.Uset2 = ones(size(biologicalData.Urxns2set));

        % Sequential gap fill
        params.sequential = 1;
        tic
        [modelList_seq] = build_network(seed_plus_iPAU,biologicalData,params);
        stseq1 = toc;

        if numel(modelList_seq{1}) > 0
            numericalIssues = 0;
        end
    end
    
    % Global gap fill
    params.sequential = 0;
    tic
    [modelList_glob] = build_network(seed_plus_iPAU,biologicalData,params);
    stseq2 = toc;
    
     % Jaccard similarity
    jaccard_sim_seq = jaccardSim(modelList_seq{1}.rxns,modPA14_v24.rxns);
    jaccard_sim_glob = jaccardSim(modelList_glob{1}.rxns,modPA14_v24.rxns);
    fprintf(['\tJaccard sim seq. = ' num2str(jaccard_sim_seq) '\t']);
    fprintf(['\tJaccard sim glob.= ' num2str(jaccard_sim_glob) '\n']);
    
    sims(i,1) = jaccard_sim_seq;
    sims(i,2) = jaccard_sim_glob;
    sims(i,3) = length(modelList_seq{1}.rxns);
    sims(i,4) = length(modelList_glob{1}.rxns);
    sims(i,5) = stseq1;
    sims(i,6) = stseq2;
end
dlmwrite('CE3_globalVsequential.tsv', sims, '\t');



