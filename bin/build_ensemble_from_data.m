function [m] = build_ensemble_from_data(numGCs,numNGCs,modelID,stochastic,rndSeed,verbose)


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