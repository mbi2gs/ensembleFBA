function [delModel] = simulateGeneDeletion(model,gene)
% Simulate a gene deletion by blocking reactions based on GPR rules
%
% Written by Matt Biggs, 2016

geneIndex = find(ismember(model.genes,gene));

% If the gene is not in the model, no need to delete anything
if isempty(geneIndex)
    delModel = model;
else
    rxns2delete = zeros(size(model.rxns));    
    for i = 1:length(model.rxns)
        curGPR = model.grRules{i};
        curGPR = strrep(curGPR,gene,'0'); % Delete gene
        curGPR = strrep(curGPR,'or','|'); % Change logical operators
        curGPR = strrep(curGPR,'and','&');
        curGPR = strrep(curGPR,'GAP_FILLED','1'); % Keep rxns not associated with a gene
        curGPR = strrep(curGPR,'Unknown','1');
        curGPR = strrep(curGPR,'Unassigned','1');
        curGPR = strrep(curGPR,'unassigned','1');
        curGPR = strrep(curGPR,'SPONTANEOUS','1');
        curGPR = regexprep(curGPR,'peg.[0-9]+','1'); % Keep all other genes
        curGPR = regexprep(curGPR,'PA14_[0-9]+','1'); % Keep all other genes
        
        if length(curGPR) > 0
            keepRxn = eval(curGPR);
        else
            keepRxn = 1;
        end
        
        if keepRxn == 0
            rxns2delete(i) = 1;
        end
    end
    
    rxns2delete = find(rxns2delete);
    
    delModel = model;
    delModel.S(:,rxns2delete) = [];
    delModel.rxns(rxns2delete) = [];
    delModel.rxnNames(rxns2delete) = [];
    delModel.rev(rxns2delete) = [];
    delModel.c(rxns2delete) = [];
    delModel.ub(rxns2delete) = [];
    delModel.lb(rxns2delete) = [];
    delModel.grRules(rxns2delete) = [];
    
end

end