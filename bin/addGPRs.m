function [updatedModelList] = addGPRs(modelList,Rxn_GPR_mapping)
% Add GPRs to each model (just copy-paste the SEED GPRs).
%
% Written by Matt Biggs, 2016

updatedModelList = cell(length(modelList),1);

for i = 1:length(modelList)
    mdl = modelList{i};
    mdl.grRules = repmat({''},size(mdl.rxns));
    
    spontRxns = {'rxn05064'}; % from curated PA14 model
    geneList = cell(0,1);
    
    for j = 1:length(mdl.rxns)
        curRxn = mdl.rxns{j};
        rxn_i = ismember(Rxn_GPR_mapping.rxns,curRxn);
        
        if sum(rxn_i) > 0
            GPR = Rxn_GPR_mapping.gprs{rxn_i};
            
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