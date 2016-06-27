function [m] = build_ensemble_from_data(numGCs,numNGCs,modelID,stochastic,rndSeed,verbose)

% Load universal reaction database
if ~exist('seed_rxns_mat','var')
    load seed_rxns_mat_phen2net
end

seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

%## Build biomass function (based on iMO1056):
% Substrates:
% ATP                       cpd00002	7933 (index)
% Glycine                   cpd00033	2573
% Histidine                 cpd00119	978
% Lysine                    cpd00039    2567
% Aspartate                 cpd00041    7223
% Glutamate                 cpd00023    7584
% Serine                    cpd00054    2202
% Threonine                 cpd00161    5266
% Asparagine                cpd00132    1349
% Glutamine                 cpd00053    3130
% Cysteine                  cpd00084    6520
% Proline                   cpd00129    6008
% Valine                    cpd00156    230
% Isoleucine                cpd00322    4296
% Leucine                   cpd00107    5622
% Methionine                cpd00060    6861
% Phenylalanine             cpd00066    6863
% Tyrosine                  cpd00069    6867
% Tryptophan                cpd00065    6864
% Cardiolipin               cpd12801    6246
% CoA                       cpd00010    2939
% CTP                       cpd00052    2204
% dATP                      cpd00115    982
% dCTP                      cpd00356    635
% dGTP                      cpd00241    7842
% dTTP                      cpd00357    634
% FAD                       cpd00015    2944
% Glycogen                  cpd00155    231
% GTP                       cpd00038    2566
% H2O                       cpd00001    7930
% NAD                       cpd00003    7932
% NADH                      cpd00004    7927
% NADP                      cpd00006    7929
% NADPH                     cpd00005    7926
% Phosphatidylethanolamine  cpd11456    6907
% peptidoglycan subunit     cpd16008    3356
% Phosphatidylglycerol      cpd11652    7958
% phosphatidylserine        cpd15657    6603
% Putrescine                cpd00118    979
% Spermidine                cpd00264    6545
% Succinyl-CoA              cpd00078    1865
% UDP-glucose               cpd00026    7581
% UTP                       cpd00062    6859
% 
% Products:
% ADP           cpd00008    7935
% H+            cpd00067    6862
% Pi            cpd00009    7934
% PPi           cpd00012    2937
biomassRxn = zeros(length(seed_rxns_mat.mets),1);
biomassRxn([7933,2573,978,2567,7223,7584,2202,5266,1349,3130, ...
            6520,6008,230,4296,5622,6861,6863,6867,6864,6246, ...
            2939,2204,982,635,7842,634,2944,231,2566,7930, ...
            7932,7927,7929,7926,6907,3356,7958,6603,979,6545, ...
            1865,7581,6859],1) = -1;
biomassRxn([7935,6862,7934,2937],1) = 1;

%## Define growth conditions (based on minimal media conditions from iMO1056)
% Cobalt        cpd00149    4926 (Index)
% Copper        cpd00058    2213
% Fe3           cpd10516    3321
% H+            cpd00067    6862
% H2O           cpd00001    7930
% Mg2           cpd00254    6753
% N2            cpd00528    6714
% NH4           cpd00013    2938
% O2            cpd00007    7928
% Pi            cpd00009    7934
% SO4           cpd00048    7222
%%% No reactions use: Cd2+          cpd01012     
%%% No reactions use: Fe2           cpd00021    
%%% No reactions use: K+            cpd00205    
%%% No reactions use: Mn2           cpd00030    
%%% No reactions use: Na            cpd00971    
%%% No reactions use: Zn2           cpd00034    
% Growth Carbon sources:
% 4-Hydroxybenzoate     cpd00136    1346
% Acetate               cpd00029    7578
% Citrate               cpd00137    1345
% Glycerol-3-phosphate  cpd00080    6524
% D-Alanine             cpd00117    980
% D-Fructose            cpd00082    6526
% Adenosine             cpd00182    4234
% Fumarate              cpd00106    5621
% L-Arginine            cpd00051    3141
% L-Asparagine          cpd00132    1349
% L-Aspartate           cpd00041    7223
% L-Glutamate           cpd00023    7584
% L-Glutamine           cpd00053    3130
% L-Histidine           cpd00119    978
% Glycerol              cpd00100    5623
% Glycine               cpd00033    2573
% Itaconate             cpd00380    8748 % Infesible
% Ornithine             cpd00064    6865
% L-Phenylalanine       cpd00066    6863
% L-Proline             cpd00129    6008
% L-Serine              cpd00054    2202
% Malonate              cpd00308    8195
% L-Alanine             cpd00035    2568
% Glucose               cpd00027    7580
% Putrescine            cpd00118    979
% Pyruvate              cpd00020    7586
% Succinate             cpd00036    2569
% 2-Oxoglutarate        cpd00024    7583
% GABA                  cpd00281    8288
% L-Isoleucine          cpd00322    4296
% L-Malate              cpd00130    1347
% L-Lactate             cpd00159    233
% L-Leucine             cpd00107    5622
% D-Mannitol            cpd00314    9684
% Aminoethanol          cpd00162    5263
% N-Acetyl-L-glutamate  cpd00477    3940
% gly-glu-L             cpd11592    3895
% gly-pro-L             cpd11588    8803
% Hydroxy-L-proline     cpd00851    1958
% Butyrate              cpd00211    1874
% L-alanylglycine       cpd11585    3450
% L-Lysine              cpd00039    2567
% Carnitine             cpd00266    6547
% GLCN                  cpd00222    7244
% Hydroxyphenylacetate  cpd00489    247
% Propionate            cpd00141    4918
% Uridine               cpd00249    6209
% Hydroxybutanoate      cpd00797    7608

lb_minimal_base = zeros(length(seed_rxns_mat.mets),1);
lb_minimal_base([4926,2213,3321,6862,7930,6753,6714,2938,7928,7934,7222],1) = -1000;

growth_carbon_sources = [1346,7578,1345,6524,980,6526,4234,5621,3141,1349, ...
                         7223,7584,3130,978,5623,2573,6865,6863,6008, ...
                         2202,8195,2568,7580,979,7586,2569,7583,8288,4296,...
                         1347,233,5622,9684,5263,3940,3895,8803,1958,1874, ...
                         3450,2567,6547,7244,247,4918,6209,7608];
               
n = length(growth_carbon_sources);
growthConditions = repmat(lb_minimal_base,[1,n]);
for i = 1:n
    growthConditions(growth_carbon_sources(i),i) = -10;
end

% Define non-growth conditions
% non-Growth Carbon sources:
% 2,3-Butanediol        cpd01949    8832
% ACTN                  cpd00361    462
% CELB                  cpd00158    234
% fructose-6-phosphate  cpd00072    1859
% D-Galacturonate       cpd00280    8289
% Glucose-1-phosphate   cpd00089    6518
% glucose-6-phosphate   cpd00079    1866
% D-Glucarate           cpd00609    6445
% D-Mannose             cpd00138    1351
% Glucuronate           cpd00164    5261
% Sorbitol              cpd00588    7827
% Tartrate              cpd00666    4395
% Xylose                cpd00154    232
% Formate               cpd00047    7229
% Glycogen              cpd00155    231
% Glycolate             cpd00139    1350
% Glyoxalate            cpd00040    7224
% gly-asp-L             cpd11589    3438
% dAMP                  cpd00294    3272
% L-Methionine          cpd00060    6861
% Thyminose             cpd01242    3171
% Thymidine             cpd00184    4229
% TRHL                  cpd00794    7607
% L-Valine              cpd00156    230
% Maltose               cpd00179    602
% Amylotriose           cpd01262    3489
% Hydroxyphenylacetate  cpd03320    519
% 2-Oxobutyrate         cpd00094    1523
% L-Inositol            cpd00121    6001
% D-Mucic acid          cpd00652    8768
% L-Arabinose           cpd00224    7238
% L-Homoserine          cpd00227    7239
% SALC                  cpd00599    3162
% Citraconate           cpd01502    3003
% Acetoacetate          cpd00142    4919
% D-Malate              cpd00386    6024
% D-Ribose              cpd00105    5620
% D-Serine              cpd00550    2431
% Inosine               cpd00246    6219
% L-Threonine           cpd00161    5266

nongrowth_carbon_sources = [8832,462,234,1859,8289,6518,1866,6445,1351,5261, ...
                            7827,4395,232,7229,231,1350,7224,3438,3272,6861, ...
                            3171,4229,7607,230,602,3489,519,1523,6001,8768, ...
                            7238,7239,3162,3003,4919,6024,5620,2431,6219,5266];

nonGrowthConditions = repmat(lb_minimal_base,[1,length(nongrowth_carbon_sources)]);
for i = 1:length(nongrowth_carbon_sources)
    nonGrowthConditions(nongrowth_carbon_sources(i),i) = -10;
end


%------------------------------------------------------------------------
% Simplify the universal rxn database so that there are no blocked rxns
%------------------------------------------------------------------------
% [consistentUnivRxnSet,keepURxns,keepXRxns,newBiomassFn] = removeBlockedReactionsAndBiomassComponents(seed_mat_small,[growthConditions nonGrowthConditions],biomassRxn);
% save('consistentUnivRxnSet.mat','consistentUnivRxnSet','keepURxns','keepXRxns','newBiomassFn')
load consistentUnivRxnSet.mat

trimmed_growthConditions = growthConditions(keepXRxns,:);
trimmed_nonGrowthConditions = nonGrowthConditions(keepXRxns,:);

%------------------------------------------------------------------------
% Import reactions from P. aeruginosa PA1 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
% fig|1279007.3.peg.2761
% rxn02201_c0
fid = fopen('PA14_reference_genome.rxntbl','r');
rxns = textscan(fid, '%s%s%s%s%s','Delimiter','\t');
fclose(fid);

trimmed_rxnList = char(rxns{1});
trimmed_rxnList = trimmed_rxnList(:,1:8);
gprs = strfind(rxns{5},'fig');
rxns_with_genomic_evidence = ~cellfun(@isempty,gprs);
trimmed_rxnList = cellstr(trimmed_rxnList(rxns_with_genomic_evidence,:));
trimmed_gprs = rxns{5};
trimmed_gprs = cellfun(@(orig,old,new) strrep(orig,old,new), trimmed_gprs, repmat({'fig|208963.12.'},[length(trimmed_gprs),1]), repmat({''},[length(trimmed_gprs),1]),'UniformOutput',false);
trimmed_gprs = trimmed_gprs(rxns_with_genomic_evidence);

Urxns2set = [find(ismember(consistentUnivRxnSet.rxns,trimmed_rxnList)); find(ismember(consistentUnivRxnSet.rxns,'rxn05064'))]; % include spontaneous rxn05064
Uset2 = ones(size(Urxns2set));

% Include exchange reactions for all non-growth conditions, just so that
% it's the network itself--not the lack of exchange reactions--that prevents growth
Xrxns2set = find(sum( abs(consistentUnivRxnSet.X([growth_carbon_sources nongrowth_carbon_sources],:)) ,1) > 0);
Xset2 = ones(size(Xrxns2set));

%------------------------------------------------------------------------
% Create a single consistent model and store in a .mat file
%------------------------------------------------------------------------
rng(rndSeed,'twister');

training_growthConditions = trimmed_growthConditions(:,1:20);
training_nonGrowthConditions = trimmed_nonGrowthConditions(:,1:20);
testing_growthConditions = trimmed_growthConditions(:,21:end);
testing_nonGrowthConditions = trimmed_nonGrowthConditions(:,21:end);

trainSetGC = [];
testSetGC = [];
if numGCs <= size(training_growthConditions,2)
    p1 = randperm(size(training_growthConditions,2));
    trainSetGC = p1(1:numGCs);
    testSetGC = 21:30; % Pick 10 GCs that weren't explicitely gap filled against as test set
    training_growthConditions = training_growthConditions(:,trainSetGC);
end

trainSetNGC = [];
testSetNGC = [];
if numNGCs <= size(training_nonGrowthConditions,2)
    p1 = randperm(size(training_nonGrowthConditions,2));
    trainSetNGC = p1(1:numNGCs);
    testSetNGC = 21:30;
    training_nonGrowthConditions = training_nonGrowthConditions(:,trainSetNGC);
end

% Random subset of gene annotations
p2 = randperm(length(Urxns2set));
len80p = ceil(length(Urxns2set) * 0.8);
Urxns2set = Urxns2set(p2(1:len80p));
Uset2 = Uset2(p2(1:len80p));

% Set parameters
biologicalData = struct;
biologicalData.growthConditions = training_growthConditions;
biologicalData.nonGrowthConditions = training_nonGrowthConditions;
biologicalData.biomassFn = newBiomassFn;
biologicalData.Urxns2set = Urxns2set;
biologicalData.Uset2 = Uset2;
biologicalData.Xrxns2set = Xrxns2set;
biologicalData.Xset2 = Xset2;

params = struct;
params.sequential = 1;
params.stochast = stochastic;
params.rndSeed = rndSeed*10;
params.numModels2gen = 1;
params.verbose = verbose;

tic
[modelList1] = phen2net_iterative_dec(consistentUnivRxnSet,biologicalData,params);
time2run = toc;

% Add GPRs
Rxn_GPR_mapping = struct;
Rxn_GPR_mapping.rxns = trimmed_rxnList;
Rxn_GPR_mapping.gprs = trimmed_gprs;
[updatedModelList] = addGPRs(modelList1,Rxn_GPR_mapping);

m = updatedModelList{1};
m.time2run = time2run;

if verbose > 0
    time2run
end

% Check which conditions the model grows in
[gc_bm_vals,ngc_bm_vals] = testModelsInGrowthConditions_flex({m},consistentUnivRxnSet.Ex_names,growthConditions(keepXRxns,:),nonGrowthConditions(keepXRxns,:),verbose);
m.gc_bm_vals = gc_bm_vals;
m.ngc_bm_vals = ngc_bm_vals;
m.trainingGCs = trainSetGC;
m.testGCs = testSetGC;
m.trainingNGCs = trainSetNGC;
m.testNGCs = testSetNGC;


save(modelID,'m');

end