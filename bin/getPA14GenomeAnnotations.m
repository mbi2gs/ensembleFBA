function [PA14GenomicData] = getPA14GenomeAnnotations()
%-------------------------------------------------------------------------- 
% getPA14GenomeAnnotations - Reads in genome annotations, peg-to-rxn
% mappings, etc.
%
% Inputs:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%       matrix for exchange rxns (X matrix), a reversability indicator for
%       all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange 
%       rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames), 
%       and metabolite formulas (metFormulas)
%
% Outputs:
%     PA14GenomicData - a Matlab struct with the following fields:
%       rxn_GPR_mapping - a Matlab struct with the following fields:
%           rxns - a cell array of rxn IDs of the same size as "rxn_GPR_mapping.gprs"
%           gprs - a cell array of GPRs of the same size as "rxn_GPR_mapping.rxns"
%
% Written by Matt Biggs
%--------------------------------------------------------------------------

%------------------------------------------------------------------------
% Import reactions from P. aeruginosa PA14 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
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

rxn_GPR_mapping = struct;
rxn_GPR_mapping.rxns = trimmed_rxnList;
rxn_GPR_mapping.gprs = trimmed_gprs;

% Export
PA14GenomicData = struct;
PA14GenomicData.rxn_GPR_mapping = rxn_GPR_mapping;

end