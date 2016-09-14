function [StrepGenomicData] = getStrepGenomeAnnotations()
%-------------------------------------------------------------------------- 
% getStrepGenomeAnnotations - Reads in genome annotations, peg-to-rxn
% mappings, etc.
%
% Inputs:
%
% Outputs:
%     StrepGenomicData - a Matlab struct with the following fields:
%        mitis_rxn_GPR_mapping, gallolyticus_rxn_GPR_mapping, oralis_rxn_GPR_mapping,
%        equinus_rxn_GPR_mapping, pneumoniae_rxn_GPR_mapping,
%        vestibularis_rxn_GPR_mapping - Matlab structs with the following fields:
%               rxns - a cell array of rxn IDs of the same size as "rxn_GPR_mapping.gprs"
%               gprs - a cell array of GPRs of the same size as "rxn_GPR_mapping.rxns"
%
% Written by Matt Biggs
%--------------------------------------------------------------------------

%------------------------------------------------------------------------
% Import reactions from Streptococcus mitis ATCC 6249 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('Streptococcus_mitis_ATCC_6249.rxntbl','r');
rxns = textscan(fid, '%s%s%s%s%s','Delimiter','\t');
fclose(fid);

trimmed_rxnList = char(rxns{1});
trimmed_rxnList = trimmed_rxnList(:,1:8);
gprs = strfind(rxns{5},'fig');
rxns_with_genomic_evidence = ~cellfun(@isempty,gprs);
trimmed_rxnList = cellstr(trimmed_rxnList(rxns_with_genomic_evidence,:));
trimmed_gprs = rxns{5};
trimmed_gprs = cellfun(@(orig,old,new) strrep(orig,old,new), trimmed_gprs, repmat({'fig|864567.3.'},[length(trimmed_gprs),1]), repmat({''},[length(trimmed_gprs),1]),'UniformOutput',false);
trimmed_gprs = trimmed_gprs(rxns_with_genomic_evidence);

rxn_GPR_mapping = struct;
rxn_GPR_mapping.rxns = trimmed_rxnList;
rxn_GPR_mapping.gprs = trimmed_gprs;

% Export
StrepGenomicData = struct;
StrepGenomicData.mitis_rxn_GPR_mapping = rxn_GPR_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus gallolyticus strain ICDDRB-NRC-S3 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('Streptococcus_gallolyticus_strain_ICDDRB-NRC-S3.rxntbl','r');
rxns = textscan(fid, '%s%s%s%s%s','Delimiter','\t');
fclose(fid);

trimmed_rxnList = char(rxns{1});
trimmed_rxnList = trimmed_rxnList(:,1:8);
gprs = strfind(rxns{5},'fig');
rxns_with_genomic_evidence = ~cellfun(@isempty,gprs);
trimmed_rxnList = cellstr(trimmed_rxnList(rxns_with_genomic_evidence,:));
trimmed_gprs = rxns{5};
trimmed_gprs = cellfun(@(orig,old,new) strrep(orig,old,new), trimmed_gprs, repmat({'fig|315405.10.'},[length(trimmed_gprs),1]), repmat({''},[length(trimmed_gprs),1]),'UniformOutput',false);
trimmed_gprs = trimmed_gprs(rxns_with_genomic_evidence);

rxn_GPR_mapping = struct;
rxn_GPR_mapping.rxns = trimmed_rxnList;
rxn_GPR_mapping.gprs = trimmed_gprs;

% Export
StrepGenomicData.gallolyticus_rxn_GPR_mapping = rxn_GPR_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus oralis ATCC 49296 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('Streptococcus_oralis_ATCC_49296.rxntbl','r');
rxns = textscan(fid, '%s%s%s%s%s','Delimiter','\t');
fclose(fid);

trimmed_rxnList = char(rxns{1});
trimmed_rxnList = trimmed_rxnList(:,1:8);
gprs = strfind(rxns{5},'fig');
rxns_with_genomic_evidence = ~cellfun(@isempty,gprs);
trimmed_rxnList = cellstr(trimmed_rxnList(rxns_with_genomic_evidence,:));
trimmed_gprs = rxns{5};
trimmed_gprs = cellfun(@(orig,old,new) strrep(orig,old,new), trimmed_gprs, repmat({'fig|888049.3.'},[length(trimmed_gprs),1]), repmat({''},[length(trimmed_gprs),1]),'UniformOutput',false);
trimmed_gprs = trimmed_gprs(rxns_with_genomic_evidence);

rxn_GPR_mapping = struct;
rxn_GPR_mapping.rxns = trimmed_rxnList;
rxn_GPR_mapping.gprs = trimmed_gprs;

% Export
StrepGenomicData.oralis_rxn_GPR_mapping = rxn_GPR_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus equinus strain AG46 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('Streptococcus_equinus_strain_AG46.rxntbl','r');
rxns = textscan(fid, '%s%s%s%s%s','Delimiter','\t');
fclose(fid);

trimmed_rxnList = char(rxns{1});
trimmed_rxnList = trimmed_rxnList(:,1:8);
gprs = strfind(rxns{5},'fig');
rxns_with_genomic_evidence = ~cellfun(@isempty,gprs);
trimmed_rxnList = cellstr(trimmed_rxnList(rxns_with_genomic_evidence,:));
trimmed_gprs = rxns{5};
trimmed_gprs = cellfun(@(orig,old,new) strrep(orig,old,new), trimmed_gprs, repmat({'fig|1335.9.'},[length(trimmed_gprs),1]), repmat({''},[length(trimmed_gprs),1]),'UniformOutput',false);
trimmed_gprs = trimmed_gprs(rxns_with_genomic_evidence);

rxn_GPR_mapping = struct;
rxn_GPR_mapping.rxns = trimmed_rxnList;
rxn_GPR_mapping.gprs = trimmed_gprs;

% Export
StrepGenomicData.equinus_rxn_GPR_mapping = rxn_GPR_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus pneumoniae genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('Streptococcus_pneumoniae.rxntbl','r');
rxns = textscan(fid, '%s%s%s%s%s','Delimiter','\t');
fclose(fid);

trimmed_rxnList = char(rxns{1});
trimmed_rxnList = trimmed_rxnList(:,1:8);
gprs = strfind(rxns{5},'fig');
rxns_with_genomic_evidence = ~cellfun(@isempty,gprs);
trimmed_rxnList = cellstr(trimmed_rxnList(rxns_with_genomic_evidence,:));
trimmed_gprs = rxns{5};
trimmed_gprs = cellfun(@(orig,old,new) strrep(orig,old,new), trimmed_gprs, repmat({'fig|1313.5731.'},[length(trimmed_gprs),1]), repmat({''},[length(trimmed_gprs),1]),'UniformOutput',false);
trimmed_gprs = trimmed_gprs(rxns_with_genomic_evidence);

rxn_GPR_mapping = struct;
rxn_GPR_mapping.rxns = trimmed_rxnList;
rxn_GPR_mapping.gprs = trimmed_gprs;

% Export
StrepGenomicData.pneumoniae_rxn_GPR_mapping = rxn_GPR_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus vestibularis strain 22-06 S6 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('Streptococcus_vestibularis_strain_22-06-S6.rxntbl','r');
rxns = textscan(fid, '%s%s%s%s%s','Delimiter','\t');
fclose(fid);

trimmed_rxnList = char(rxns{1});
trimmed_rxnList = trimmed_rxnList(:,1:8);
gprs = strfind(rxns{5},'fig');
rxns_with_genomic_evidence = ~cellfun(@isempty,gprs);
trimmed_rxnList = cellstr(trimmed_rxnList(rxns_with_genomic_evidence,:));
trimmed_gprs = rxns{5};
trimmed_gprs = cellfun(@(orig,old,new) strrep(orig,old,new), trimmed_gprs, repmat({'fig|1343.5.'},[length(trimmed_gprs),1]), repmat({''},[length(trimmed_gprs),1]),'UniformOutput',false);
trimmed_gprs = trimmed_gprs(rxns_with_genomic_evidence);

rxn_GPR_mapping = struct;
rxn_GPR_mapping.rxns = trimmed_rxnList;
rxn_GPR_mapping.gprs = trimmed_gprs;

% Export
StrepGenomicData.vestibularis_rxn_GPR_mapping = rxn_GPR_mapping;

end