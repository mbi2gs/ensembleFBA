function [ensemble] = build_ensemble(universalRxnSet,biologicalData,params)
%-------------------------------------------------------------------------- 
% build_ensemble - Generates a set of metabolic network reconstructions
% trained on subsets of the data and/or permutations of the growth
% conditions.
%
% Inputs:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%       matrix for exchange rxns (X matrix), a reversability indicator for
%       all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange 
%       rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames), 
%       and metabolite formulas (metFormulas)
%     biologicalData - Matlab structure with several optional elements:
%           growthConditions = set of lower bounds corresponding to growth media conditions
%           nonGrowthConditions = set of lower bounds corresponding to non-growth media conditions
%           biomassFn = same format as a reaction in universalRxnSet.S
%           Urxns2set = A list of rxn indices (in universalRxnSet.S) that are forced to be included or excluded
%           Uset2 = List of 1's or 0's (inclusion or exclusion)
%           Xrxns2set = Same as Urxns2set, but for exchange rxns (universalRxnSet.X)
%           Xset2 = Same as Uset2 but for exchange rxns (universalRxnSet.X)
%           rxn_GPR_mapping = Matlab structure containing the fields:
%               rxns = a cell array of reaction IDs (as strings)
%               gprs = a cell array of gene-protein-reaction relationships (as strings)
%     params - Matlab structure with several optional elements:
%           sequential = 1 (default) indicates sequential gap filling, 0 indicates global
%           rndSequence = 1 (default) indicates a randomized order for sequential gap filling
%           fractionUrxns2set = The fraction of rxns from Urxns2set which will be used to reconstruct each network (default 0.8)
%           stochast = 1 (default) indicates stochastic weights during the expansion step, 0 indicates equal weights
%           rndSeed = Allows the user to manually set the random seed (is set to 0 as default)
%           numModels2gen = Indicates the number of models to produce (1 is the default)
%           verbose = 1 (default) indicates verbose output, 0 indicates that messages during runtime are not displayed
%
% Outputs:
%     modelList - a cell array of COBRA-format models (Matlab structs)
%
% Written by Matt Biggs
%--------------------------------------------------------------------------

% Unpack parameters
if isfield(params,'verbose')
    verbose = params.verbose;
else
    verbose = 0;
end

if isfield(params,'numModels2gen')
    numModels2gen = params.numModels2gen;
    if numModels2gen <= 0
       numModels2gen = 1; 
    end
else
    numModels2gen = 1;
end

if isfield(params,'stochast')
    stochastic = params.stochast;
else
    stochastic = 0;
end

if isfield(params,'rndSequence')
    rndSequence = params.rndSequence;
else
    rndSequence = 1;
end

if isfield(params,'fractionUrxns2set')
    fractionUrxns2set = params.fractionUrxns2set;
    
    if fractionUrxns2set < 0 || fractionUrxns2set > 1
        error('The parameter "fractionUrxns2set" is outside the range (0,1].');
    end
else
    fractionUrxns2set = 0.8;
end

if isfield(params,'rndSeed')
    rndSeed = params.rndSeed;
else
    rndSeed = now;
end

if isfield(params,'sequential')
    sequential = params.sequential;
else
    sequential = 1;
end

% Make sure biological data is present
if isfield(biologicalData,'growthConditions')
    growthConditions = biologicalData.growthConditions;
else
    error('No growth conditions were provided.');
end

if isfield(biologicalData,'nonGrowthConditions')
    nonGrowthConditions = biologicalData.nonGrowthConditions;
else
    error('No non-growth conditions were provided.');
end

if isfield(biologicalData,'biomassFn')
    biomassFn = biologicalData.biomassFn;
else
    error('No biomass objective was provided.');
end

if isfield(biologicalData,'Urxns2set')
    Urxns2set = biologicalData.Urxns2set;
    
    if isfield(biologicalData,'Uset2')
        Uset2 = biologicalData.Uset2;
    else
        error('No reaction state was provided for "Urxns2set".');
    end
    
    if length(find(Urxns2set)) ~= length(find(Uset2))
        error('"Urxns2set" and "Uset2" are different lengths.');
    end
else
    Urxns2set = [];
    Uset2 = [];
end

if isfield(biologicalData,'Xrxns2set')
    Xrxns2set = biologicalData.Xrxns2set;
    
    if isfield(biologicalData,'Xset2')
        Xset2 = biologicalData.Xset2;
    else
        error('No reaction state was provided for "Xrxns2set".');
    end
    
    if length(find(Xrxns2set)) ~= length(find(Xset2))
        error('"Xrxns2set" and "Xset2" are different lengths.');
    end
else
    Xrxns2set = [];
    Xset2 = [];
end

if isfield(biologicalData,'rxn_GPR_mapping')
    rxn_GPR_mapping = biologicalData.rxn_GPR_mapping;
    
    if ~isfield(rxn_GPR_mapping,'rxns')
        error('Missing reaction list in "rxn_GPR_mapping.rxns".');
    end
    
    if ~isfield(rxn_GPR_mapping,'gprs')
        error('Missing GPR list in "rxn_GPR_mapping.gprs".');
    end
    
    if length(rxn_GPR_mapping.rxns) ~= length(rxn_GPR_mapping.gprs)
        error('"rxn_GPR_mapping.rxns" and "rxn_GPR_mapping.gprs" are different lengths.');
    end
else
    rxn_GPR_mapping = struct;
end

%------------------------------------------------------------------------
% Create a set of consistent networks and store as a cell array of Matlab
% structures
%------------------------------------------------------------------------
rng(rndSeed,'twister');
ensemble = cell(numModels2gen,1);

if verbose > 0
    fprintf('Starting to build ensemble.\n');
end

for i = 1:numModels2gen
    % Generate Random subset of gene annotations
    if fractionUrxns2set < 1
        p1 = randperm(length(Urxns2set));
        len80p = ceil(length(Urxns2set) * fractionUrxns2set);
        tmp_Urxns2set = Urxns2set(p1(1:len80p));
        tmp_Uset2 = Uset2(p1(1:len80p));
    else
        tmp_Urxns2set = Urxns2set;
        tmp_Uset2 = Uset2;
    end
    
    % Generate random permuation of growth conditions
    if rndSequence > 0
        p1 = randperm(size(growthConditions,2));
        tmp_growthConditions = growthConditions(:,p1);
        
        p2 = randperm(size(nonGrowthConditions,2));
        tmp_nonGrowthConditions = nonGrowthConditions(:,p2);
    else
        tmp_growthConditions = growthConditions;
        tmp_nonGrowthConditions = nonGrowthConditions;
    end
    
    % Set parameters
    biologicalData_inner = struct;
    biologicalData_inner.growthConditions = tmp_growthConditions;
    biologicalData_inner.nonGrowthConditions = tmp_nonGrowthConditions;
    biologicalData_inner.biomassFn = biomassFn;
    biologicalData_inner.Urxns2set = tmp_Urxns2set;
    biologicalData_inner.Uset2 = tmp_Uset2;
    biologicalData_inner.Xrxns2set = Xrxns2set;
    biologicalData_inner.Xset2 = Xset2;
    
    params_inner = struct;
    params_inner.sequential = sequential;
    params_inner.stochast = stochastic;
    params_inner.rndSeed = rndSeed;
    params_inner.numModels2gen = 1;
    params_inner.verbose = verbose;
    
    % Build network
    tic
    [modelList] = build_network(universalRxnSet,biologicalData_inner,params_inner);
    time2run = toc;
    
    % Add GPRs
    if isfield(rxn_GPR_mapping,'rxns')
        [modelList] = addGPRs(modelList,rxn_GPR_mapping);
    end
    
    m = modelList{1};
    m.time2run = time2run;
    
    ensemble{i} = m;
    
    if verbose > 0
        fprintf('Network number %d took %d seconds to build.\n',i,time2run);
    end
end

if verbose > 0
    fprintf('Completed building ensemble.\n');
end

end