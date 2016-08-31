% Computational Experiment
% Explore how individual GENREs in an ensemble differ from each other in
% terms of coverage of the gold-standard iPAU1129
%
% Written by Matt Biggs, 2016

% Load universal reaction database and add exchange rxns
load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Get the PA14 data formatted to work with the SEED database
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
params.stochast = 1;
params.numModels2gen = 1;
params.verbose = 0;
params.sequential = 1;

%------------------------------------------------------------------------
% Gap fill sequentially, in different orders
%------------------------------------------------------------------------
N_iter = 100;
N_gcs = 25;
N_rxns = ceil(length(iPAU_nonExchangeRxns)*0.8);

initialRxnDistribution = cell(N_iter,1);
ensembleS = cell(N_iter,1);
for i = 1:N_iter
    fprintf(['\titeration ' num2str(i) ':\n']);

    numericalIssues = 1;
    while numericalIssues == 1
        c=clock; 
        rng(c(6));
        % Force networks to contain a subset of PA14 reactions
        rxnSubset = randperm(length(iPAU_nonExchangeRxns),N_rxns);
        tmpRxnList = iPAU_nonExchangeRxns(rxnSubset);
        biologicalData.Urxns2set = find(ismember(seed_plus_iPAU.rxns,tmpRxnList));
        biologicalData.Uset2 = ones(size(biologicalData.Urxns2set));
        
        % Randomly select growth conditions
        rp = randperm(size(growthConditions,2),N_gcs);
        biologicalData.growthConditions = growthConditions(:,rp);
        
        % Stochastic        
        tic
        [modelList_seq] = build_network(seed_plus_iPAU,biologicalData,params);
        stseq1 = toc;
        
        ensembleS{i} = modelList_seq{1};
        initialRxnDistribution{i} = tmpRxnList;
        
        if numel(modelList_seq{1}) > 0
            numericalIssues = 0;
        end
    end
end
save('CE12_ensembles.mat','ensembleS');

load CE12_ensembles.mat

%----------------------------------------------------------
% Check the distribution of reactions throughout ensembleS
%----------------------------------------------------------
u_rxnsS = {};
for i = 1:N_iter
    u_rxnsS = [u_rxnsS; ensembleS{i}.rxns];
end
u_rxnsS = [u_rxnsS; iPAU_nonExchangeRxns];
u_rxnsS = setdiff(u_rxnsS,seed_plus_iPAU.Ex_names);
u_rxnsS = unique(u_rxnsS);

% Count the number of times each "incorrect" rxn occurs in ensemble
noniPAUrxnsS = setdiff(u_rxnsS,iPAU_nonExchangeRxns);
countsS = zeros(size(noniPAUrxnsS));
for i = 1:length(noniPAUrxnsS)
    for j = 1:N_iter
       im = sum(ismember(ensembleS{j}.rxns,noniPAUrxnsS{i}));
       if im > 0
          countsS(i) = countsS(i) + 1; 
       end
    end
end

% Count the number of times each "correct" rxn occurs in ensemble
counts2S = zeros(size(iPAU_nonExchangeRxns));
for i = 1:length(iPAU_nonExchangeRxns)
    for j = 1:N_iter
       im = sum(ismember(ensembleS{j}.rxns,iPAU_nonExchangeRxns{i}));
       if im > 0
          counts2S(i) = counts2S(i) + 1; 
       end
    end
end

% Check the distribution of "correct" rxns initially
counts3S = zeros(size(iPAU_nonExchangeRxns));
for i = 1:length(iPAU_nonExchangeRxns)
    for j = 1:N_iter
       im = sum(ismember(initialRxnDistribution{j},iPAU_nonExchangeRxns{i}));
       if im > 0
          counts3S(i) = counts3S(i) + 1; 
       end
    end
end

% Plot the distribution of "correct" rxns together with "incorrect"
% Where correct means "found in iPAU1129"
barPlotDataS = [histc(counts2S,0:N_iter) histc(countsS,0:N_iter) histc(counts3S,0:N_iter)];
figure(1); bar(0:N_iter,barPlotDataS(:,1:2),'stacked');
figure(2);bar(0:N_iter,barPlotDataS(:,3),'stacked');
dlmwrite('CE12_stochastic_rxn_dist.tsv',barPlotDataS,'\t');

% Plot proportion of correct rxns by frequency
proportions = barPlotDataS(2:end,1) ./ sum(barPlotDataS(2:end,1:2),2);
figure(3); bar(1:N_iter,proportions,'stacked');

% Write the list of rxns that are inferred correctly in 100% of GENREs
CalwaysInferred = iPAU_nonExchangeRxns(counts2S == 100);
fid = fopen('CE12_correct_rxns_always_inferred.tsv','w');
for i = 1:length(CalwaysInferred)
    fprintf(fid,[CalwaysInferred{i} '\n']);
end
fclose(fid);

IalwaysInferred = noniPAUrxnsS(countsS == 100);
fid = fopen('CE12_incorrect_rxns_always_inferred.tsv','w');
for i = 1:length(IalwaysInferred)
    fprintf(fid,[IalwaysInferred{i} '\n']);
end
fclose(fid);




