function [solutions] = ensembleFBA(ensemble,exchangeRxnsIDs,conditions,verbose)
%-------------------------------------------------------------------------- 
% ensembleFBA - Solves a flux balanc analysis (FBA) problem for each member
% of the ensemble and each set of lower bounds.
%
% Inputs:
%     ensemble - Cell array of metabolic network reconstructions in COBRA
%       format
%     exchangeRxnsIDs - A cell array of exchange reaction IDs which
%       correspond to the rows of the "conditions" matrix
%     conditions - Each column is a set of lower bounds for the reactions
%       listed in "exchangeRxnsIDs"
%     verbose - 1 indicates that information about the FBA function will be
%       shown
%
% Outputs:
%     solutions - A matrix of size (NumColumns in Conditions x NumMembersInEnsemble).
%       Each entry indicates the flux through biomass for the model in
%       column i and the growth condition in row j.
%
% Written by Matt Biggs
%--------------------------------------------------------------------------

solutions = zeros(size(conditions,2),length(ensemble));

for i = 1:length(ensemble)
    for j = 1:size(conditions,2)
        [growth,~] = fba_flex(ensemble{i},exchangeRxnsIDs,conditions,verbose);
        solutions(j,i) = growth;
    end
end

end