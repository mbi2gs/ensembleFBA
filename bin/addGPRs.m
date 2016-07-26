function [updatedModelList] = addGPRs(modelList,rxn_gpr_mapping)
%-------------------------------------------------------------------------- 
% addGPRs - Add gene-protein-reaction associations to each model in the
% list.
%
% Inputs:
%     modelList - a cell array of COBRA-format models (Matlab structs)
%     rxn_GPR_mapping - a Matlab struct with the following fields:
%           rxns - a cell array of rxn IDs of the same size as "rxn_GPR_mapping.gprs"
%           gprs - a cell array of GPRs of the same size as "rxn_GPR_mapping.rxns"
%
% Outputs:
%     updatedModelList - a cell array of COBRA-format models (Matlab structs)
%           now with genes and GPR information
%
% Written by Matt Biggs, mb3ad@virginia.edu, 2016
%-------------------------------------------------------------------------- 

updatedModelList = cell(length(modelList),1);

for i = 1:length(modelList)
    mdl = modelList{i};
    mdl.grRules = repmat({''},size(mdl.rxns));
    
    spontRxns = {'rxn05064'}; % from curated PA14 model
    geneList = cell(0,1);
    
    for j = 1:length(mdl.rxns)
        curRxn = mdl.rxns{j};
        rxn_i = ismember(rxn_gpr_mapping.rxns,curRxn);
        
        if sum(rxn_i) > 0
            GPR = rxn_gpr_mapping.gprs{rxn_i};
            
            % Extract all gene IDs
            gs = regexp(GPR,'peg.\d+','match');
            geneList = [geneList gs];
        else
            GPR = 'GAP_FILLED';
        end
        
        rxn_i = ismember(spontRxns,curRxn);
        if sum(rxn_i) > 0
            GPR = 'SPONTANEOUS';
        end
        
        mdl.grRules{j} = GPR;
        mdl.genes = unique(geneList(:));
    end
    
    updatedModelList{i} = mdl;
end

end