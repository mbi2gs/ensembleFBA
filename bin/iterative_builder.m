function [modelList] = iterative_builder(universalRxnSet,biologicalData,params)
%-------------------------------------------------------------------------- 
% iterative_builder - Iteratively gap fills a model by first "expanding"
% (adding reactions) so that it can produce biomass in all the growth 
% conditions, then by "trimming" (removing reactions) so that it does not
% grow in the non-growth conditions.
%
% Inputs:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%     matrix for exchange rxns (X matrix), a reversability indicator for
%     all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange 
%     rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames), 
%     and metabolite formulas (metFormulas).
%     biologicalData - Matlab structure with several optional elements:
%           growthConditions = set of lower bounds corresponding to growth media conditions
%           nonGrowthConditions = set of lower bounds corresponding to non-growth media conditions
%           biomassFn = same format as a reaction in universalRxnSet.S
%           Urxns2set = A list of rxn indices (in universalRxnSet.S) that are forced to be included or excluded
%           Uset2 = List of 1's or 0's (inclusion or exclusion)
%           Xrxns2set = Same as Urxns2set, but for exchange rxns (universalRxnSet.X)
%           Xset2 = Same as Uset2 but for exchange rxns (universalRxnSet.X)
%     params - Matlab structure with several optional elements:
%           sequential = 1 (default) indicates sequential gap filling, 0 indicates global
%           stochast = 1 (default) indicates stochastic weights on the reactions in universalRxnSet.S, 0 indicates deterministic weights
%           rndSeed = Allows the user to manually set the random seed (is set to 0 as default)
%           numModels2gen = Indicates the number of models to produce (1 is the default)
%           verbose = 1 (default) indicates verbose output, 0 indicates that messages during runtime are not displayed
%
% Outputs:
%     modelList - a cell array of COBRA-format models (Matlab structs)
%
% Written by Matt Biggs, mb3ad@virginia.edu, 2016
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
    stochast = params.stochast;
else
    stochast = 0;
end

if isfield(params,'rndSeed')
    rndSeed = params.rndSeed;
else
    rndSeed = 0;
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

jaccardSim = @(a,b) sum(ismember(a,b))/length(a)';

%----------------------------------------------------
% Run the algorithm 'numModels2gen' times
%----------------------------------------------------
modelList = cell(numModels2gen,1);
for i = 1:numModels2gen
    if verbose > 0
        fprintf(['Iteration ' num2str(i) '\n']);
    end   

    iterate = 1;
    currRxnSet = cell(0,1);
    tmpNonGrowthConditions = nonGrowthConditions;
    
    % These store a master copy of the user's rxn set list with rxns 
    % removed that were exlcuded by the trim step
    cpUrxns2set = Urxns2set;
    cpUset2 = Uset2;
    cpXrxns2set = Xrxns2set;
    cpXset2 = Xset2;
    
    % These store a working copy of the rxn set list 
    tmpUrxns2set = Urxns2set;
    tmpUset2 = Uset2;
    tmpXrxns2set = Xrxns2set;
    tmpXset2 = Xset2;

    % These store the rxns which are excluded by the trim step
    tmpUexclude = [];

    % These store the most recent changes
    lastRxnDatabase = universalRxnSet;
    lastRxnDatabase.growthConditions = growthConditions;
    lastRxnDatabase.nonGrowthConditions = nonGrowthConditions;
    lastUexclude = [];
    
    % Keep track of the most consistent model and return it
    mostConsistentMdl = struct;
    consistWithNGCs = size(nonGrowthConditions,2);
    handicap = 0;
    
    while iterate > 0
 
        %----------------------------------------------------
        % Expand
        %----------------------------------------------------
        if sequential > 0
            
            for j = 1:size(growthConditions,2)

                [mdl, rxnDatabase, ~, ~, feasible] = expand(universalRxnSet,...
                                                   growthConditions(:,j),...
                                                   tmpNonGrowthConditions,...
                                                   biomassFn,...
                                                   tmpUrxns2set,tmpUset2,...
                                                   tmpXrxns2set,tmpXset2,...
                                                   verbose,stochast,rndSeed);
                
                if feasible > 0
                    % Update working rxn set list
                    expandedUrxns2set = find( ismember(universalRxnSet.rxns, rxnDatabase.rxns) );
                    expandedUset2 = ones(size(expandedUrxns2set));
                    expandedXrxns2set = find( ismember(universalRxnSet.Ex_names, rxnDatabase.Ex_names) );
                    expandedXset2 = ones(size(expandedXrxns2set));

                    tmpUrxns2set = [expandedUrxns2set(:); cpUrxns2set(:); tmpUexclude];
                    tmpUset2 = [expandedUset2(:); cpUset2(:); zeros(size(tmpUexclude))];
                    tmpXrxns2set = [expandedXrxns2set(:); cpXrxns2set(:)];
                    tmpXset2 = [expandedXset2(:); cpXset2(:)];
                else
                    break;
                end
                
            end
            
            if feasible > 0
                rxnDatabase.growthConditions = growthConditions(ismember(universalRxnSet.Ex_names,rxnDatabase.Ex_names),:);
            end
            
        else
            [mdl, rxnDatabase, ~, ~, feasible] = expand(universalRxnSet,...
                                               growthConditions,...
                                               tmpNonGrowthConditions,...
                                               biomassFn,...
                                               tmpUrxns2set,tmpUset2,...
                                               tmpXrxns2set,tmpXset2,...
                                               verbose,stochast,rndSeed);
        end
        
        %----------------------------------------------------
        % Trim
        %----------------------------------------------------    
        if feasible > 0 
            % Find correspondence between exchange reaction indices
            ex_i_m = find(ismember(mdl.rxns,universalRxnSet.Ex_names));
            ex_i_gcs = zeros(size(ex_i_m));
            
            for j = 1:length(ex_i_m)
                curExRxn = mdl.rxns{ex_i_m(j)};
                ex_i_gcs(j) = find(ismember(universalRxnSet.Ex_names,curExRxn));
            end
            
            % Identify the reactions which are used in all growth
            % conditions
            growthRxnsList = cell(size(growthConditions,2),1);
            for j = 1:size(growthConditions,2)
                mdl.lb = zeros(size(mdl.lb));
                mdl.lb(mdl.rev > 0) = -1000;
                mdl.lb(ex_i_m) = growthConditions(ex_i_gcs,j);
                mdl.ub = 1000*ones(size(mdl.lb));
                
                % FBA
                [sol_p,~,~] = fba_phen2net(mdl,0);
                growth = sol_p.x(mdl.c > 0);
                
                rxnsCarryingFlux = find(sol_p.x);
                GCRxns = find(ismember(rxnDatabase.rxns,mdl.rxns(rxnsCarryingFlux)));
                uGCRxns = unique(GCRxns);
                
                growthRxnsList{j} = uGCRxns;
            end
            
            % Identify the reactions which are used in any non-growth
            % conditions
            sequentialTrimmedRxns = [];
            didItGrowInNGC = 0;
            for j = 1:size(tmpNonGrowthConditions,2)
                mdl.lb = zeros(size(mdl.lb));
                mdl.lb(mdl.rev > 0) = -1000;
                mdl.lb(ex_i_m) = tmpNonGrowthConditions(ex_i_gcs,j);
                mdl.ub = 1000*ones(size(mdl.lb));
                
                % FBA
                [sol_p,~,~] = fba_phen2net(mdl,0);
                growth = sol_p.x(mdl.c > 0);
                
                % If there's growth, find rxns which are carrying flux
                % Don't trim any more if there are already 10 reactions
                % being removed
                if growth > 0
                    didItGrowInNGC = didItGrowInNGC + 1;
                end
                if growth > 0 && sum(~sequentialTrimmedRxns) < 10
                    rxnsCarryingFlux = find(sol_p.x);
                    NGCRxns = find(ismember(rxnDatabase.rxns,mdl.rxns(rxnsCarryingFlux)));
                    uNGCRxns = unique(NGCRxns);
                    
                    % Find growth conditions with similar reaction usage
                    rxnUsageSim = cellfun(@(x) jaccardSim(uNGCRxns,x), growthRxnsList);
                    [~,I] = sort(rxnUsageSim,'descend');
                    
                    mostSimGCs = 1:size(growthConditions,2);
                    if size(growthConditions,2) > 4
                        mostSimGCs = I(1:5);
                    end
                    
                    % Trim reactions, and if necessary, growth conditions
                    [mdl2, ~, fromDecthresh, fromGCthresh] = trim_active(rxnDatabase, ...
                                                    rxnDatabase.growthConditions(:,mostSimGCs),...
                                                    rxnDatabase.nonGrowthConditions(:,j),...
                                                    biomassFn,...
                                                    uNGCRxns,...
                                                    verbose);
                    if isempty(sequentialTrimmedRxns)
                        sequentialTrimmedRxns = fromDecthresh;
                    elseif ~isempty(fromDecthresh)
                        sequentialTrimmedRxns = sequentialTrimmedRxns & fromDecthresh;
                    end
                end
            end
            
            % Keep track of most consistent model
            if (didItGrowInNGC + handicap) < consistWithNGCs
                mostConsistentMdl = mdl;
                consistWithNGCs = didItGrowInNGC;
            end
            
            % Update the excluded rxns
            lastRxnDatabase = rxnDatabase;
            lastUexclude = find( ismember(universalRxnSet.rxns, rxnDatabase.rxns(~sequentialTrimmedRxns)) );
            tmpUexclude = [tmpUexclude; lastUexclude(:)];
            
            % Remove conflicts from master rxn set list
            cpUr2s_i = ismember(cpUrxns2set,tmpUexclude);
            cpUrxns2set(cpUr2s_i) = [];
            cpUset2(cpUr2s_i) = [];
            
            % Update working rxn set list
            tmpUrxns2set = [cpUrxns2set(:); tmpUexclude];
            tmpUset2 = [cpUset2(:); zeros(size(tmpUexclude))];
            tmpXrxns2set = Xrxns2set;
            tmpXset2 = Xset2;
            
            % Check if model is still shrinking
            if didItGrowInNGC < 1
                iterate = 0;
            elseif (jaccardSim(mdl.rxns,currRxnSet) + jaccardSim(currRxnSet,mdl.rxns)) < 2 ...
               || isempty(currRxnSet)
                currRxnSet = mdl.rxns;
            else
                iterate = 0;
            end
            
        elseif size(tmpNonGrowthConditions,2) > 0
            % Revert rxn exclusion list to last step before infeasibility
            tmpUexclude(ismember(tmpUexclude,lastUexclude)) = [];
            
            % Revert the master rxn set list to last step before infeasibility
            cpUrxns2set = Urxns2set(~ismember(Urxns2set,tmpUexclude));
            cpUset2 = Uset2(~ismember(Urxns2set,tmpUexclude));
            
            % Update working rxn set list
            tmpUrxns2set = [cpUrxns2set(:); tmpUexclude];
            tmpUset2 = [cpUset2(:); zeros(size(tmpUexclude))];
            tmpXrxns2set = Xrxns2set;
            tmpXset2 = Xset2;
            
            % Trim inconsistent non-growth conditions
            lastRxnDatabase.rxns
            [mdl2, ~, fromNGCthresh] = trim_ngc(lastRxnDatabase, ...
                                                    lastRxnDatabase.growthConditions,...
                                                    lastRxnDatabase.nonGrowthConditions,...
                                                    biomassFn,...
                                                    [],[],...
                                                    [],[],...
                                                    verbose);
            
            tmpNonGrowthConditions = tmpNonGrowthConditions(:,fromNGCthresh);  
            
            handicap = handicap + sum(~fromNGCthresh);
            
            if verbose > 0
                numRemoved = length(fromNGCthresh) - sum(fromNGCthresh);
                fprintf(['\nRemoved ' num2str(numRemoved) ' inconsistent growth condition(s).\n\n'])
            end
            
        else
            mdl = struct;
            iterate = 0;
            fprintf('\nWARNING: There is an inconsistent growth condition.\n\n')
        end
    end

    modelList{i,1} = mostConsistentMdl;    
end

end