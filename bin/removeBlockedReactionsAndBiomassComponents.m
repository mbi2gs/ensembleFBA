function [consistentUnivRxnSet,keepURxns,keepXRxns,newBiomassFn] = removeBlockedReactionsAndBiomassComponents(universalRxnSet,allGrowthConditions,biomassFn)
% Given a "universal" reaction database (such as from Model SEED), a set 
% of growth conditions (lower bounds on exchange fluxes), and a list of
% biomass components, this script removes all reactions from the database
% which cannot carry flux under any growth condition, and removes biomass
% components which cannot be synthesized under all growth condition.
%
% Written by Matt Biggs, 2016

n_agc = size(allGrowthConditions,2);
consistentRxnsAllConditions = zeros(length(universalRxnSet.rxns)+length(universalRxnSet.mets),1);

% Get the list of reactions that can carry flux under at least one
% condition
fprintf('start finding blocked reactions\n')
tmpM = struct;
tmpM.S = [universalRxnSet.S universalRxnSet.X];
for i = 1:n_agc
    tmpM.lb = double([-1000*universalRxnSet.rev(:); allGrowthConditions(:,i)]);

    [~,consistentRxns] = findBlockedReactions(tmpM);
    
    consistentRxnsAllConditions = consistentRxnsAllConditions + consistentRxns;
end
fprintf('start reformating UX\n')

% Reformat universal rxn set by removing reactions that cannot carry flux
keepURxns = consistentRxnsAllConditions > 0;
keepURxns(size(universalRxnSet.S,2)+1:end) = [];
keepXRxns = consistentRxnsAllConditions > 0;
keepXRxns(1:size(universalRxnSet.S,2)) = [];
if sum(consistentRxnsAllConditions > 0) == length(consistentRxnsAllConditions)
    consistentUnivRxnSet = universalRxnSet;
else
    consistentUnivRxnSet = struct;
    consistentUnivRxnSet.mets = universalRxnSet.mets;
    consistentUnivRxnSet.metNames = universalRxnSet.metNames;
    consistentUnivRxnSet.metFormulas = universalRxnSet.metFormulas;
    consistentUnivRxnSet.rxns = universalRxnSet.rxns(keepURxns);
    consistentUnivRxnSet.rxnNames = universalRxnSet.rxnNames(keepURxns);
    consistentUnivRxnSet.rev = universalRxnSet.rev(keepURxns);
    consistentUnivRxnSet.S = universalRxnSet.S(:,keepURxns);
    consistentUnivRxnSet.X = universalRxnSet.X(:,keepXRxns);
    consistentUnivRxnSet.Ex_names = universalRxnSet.Ex_names(keepXRxns);
end
fprintf('done finding blocked reactions\n')
fprintf('start identifying blocked biomass components\n')

% Get the list of biomass components that can be synthesized under all
% conditions
biomassComponents2Make = find(biomassFn < 0);
bc2m_allConditions = zeros(size(biomassComponents2Make));
n_bc2m = length(biomassComponents2Make);

tmpM = struct;
tmpM.S = [universalRxnSet.S universalRxnSet.X];
tmpM.ub = 1000*ones(size(tmpM.S,2),1);
tmpM.b = zeros(size(tmpM.S,1),1);
for i = 1:n_bc2m
    tmpM.c = zeros(size(tmpM.S,2),1);
    tmpM.c(size(universalRxnSet.S,2) + biomassComponents2Make(i)) = 1;
    fprintf('biomass component\n')
    j = 1;
    canBeMade = 1;    
    while canBeMade > 0
        tmpM.lb = double([-1000*universalRxnSet.rev(:); allGrowthConditions(:,j)]);
        fprintf('condition ')
        [p,~,~] = fba_phen2net(tmpM,0);

        if p.objval == 0 || j == n_agc
            canBeMade = 0;
        end

        j = j + 1;
    end
    fprintf('\n')
    if j >= n_agc
        bc2m_allConditions(i) = 1;
    end
end
fprintf('reformatting biomass\n')

% Remove biomass components which cannot be synthesized in at least one condition
if sum(bc2m_allConditions) == length(bc2m_allConditions)
    newBiomassFn = biomassFn;
else
    newBiomassFn = zeros(size(biomassFn));
    newBiomassFn( biomassComponents2Make(bc2m_allConditions > 0) ) = biomassFn( biomassComponents2Make(bc2m_allConditions > 0) );
end
fprintf('done\n')
end