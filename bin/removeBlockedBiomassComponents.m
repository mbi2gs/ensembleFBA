function [newBiomassFn] = removeBlockedBiomassComponents(universalRxnSet,growthConditions,biomassFn)
% removeBlockedBiomassComponents - Given a "universal" reaction database 
% (such as from Model SEED), a set of growth conditions (lower bounds on 
% exchange fluxes), and a list of biomass components, this script removes 
% all biomass components which cannot be synthesized under every growth 
% condition.
%
% Inputs:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%     matrix for exchange rxns (X matrix), and a reversability indicator for
%     all rxns in S.
%     growthConditions - matrix of lower bounds with rows corresponding to
%     reactions in the universalRxnSet.X matrix and columns corresponding to
%     growth conditions under which biomass should be produced.
%     biomassFn - A vector with elements corresponding to the rows in 
%     universalRxnSet.S, where each entry is the stoichiometric coefficient for
%     a metabolite in the biomass rxn.
%
% Outputs:
%     newbiomassFn = A vector of the same format as the input "biomassFn" but
%     with potentially altered entries.
% 
% Written by Matt Biggs, 2016

n_agc = size(growthConditions,2);
fprintf('start identifying blocked biomass components\n')

% Get the list of biomass components that can be synthesized under all
% conditions
biomassComponents2Make = find(biomassFn < 0);
bc2m_allGrowthConditions = zeros(size(biomassComponents2Make));
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
        tmpM.lb = double([-1000*universalRxnSet.rev(:); growthConditions(:,j)]);
        fprintf('condition ')
        [p,~,~] = fba_phen2net(tmpM,0);
        
        if p.objval < 1e-10 || j == n_agc
            canBeMade = 0;
        end
        
        j = j + 1;
    end
    fprintf('\n')
    if j > n_agc
        bc2m_allGrowthConditions(i) = 1;
    else
        biomassComponents2Make(i)
    end
end
fprintf('reformatting biomass\n')

% Remove biomass components which cannot be synthesized in at least one condition
if sum(bc2m_allGrowthConditions) == length(bc2m_allGrowthConditions)
    newBiomassFn = biomassFn;
else
    newBiomassFn = zeros(size(biomassFn));
    newBiomassFn( biomassComponents2Make(bc2m_allGrowthConditions > 0) ) = biomassFn( biomassComponents2Make(bc2m_allGrowthConditions > 0) );
    newBiomassFn(biomassFn > 0) = biomassFn(biomassFn > 0);
end

fprintf('done\n')
end