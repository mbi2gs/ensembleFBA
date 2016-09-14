function [StrepGeneDrugPairs] = getStrepGeneDrugPairs()
%-------------------------------------------------------------------------- 
% getStrepGeneDrugPairs - Reads in gene-drug pairs as output by
% "listPossibleTargets.py".
%
% Inputs:
%
% Outputs:
%     StrepGeneDrugPairs - a Matlab struct with the following fields:
%        mitis_gene_drug_mapping, gallolyticus_gene_drug_mapping, oralis_gene_drug_mapping,
%        equinus_gene_drug_mapping, pneumoniae_gene_drug_mapping,
%        vestibularis_gene_drug_mapping - Matlab structs with the following fields:
%               genes - a cell array of gene IDs of the same size as "gene_drug_mapping.drugIDs"
%               drugIDs - a cell array of GPRs of the same size as "gene_drug_mapping.genes"
%
% Written by Matt Biggs, 2016
%--------------------------------------------------------------------------

%------------------------------------------------------------------------
% Import reactions from Streptococcus mitis ATCC 6249 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('s.mitis.results.targets','r');
geneDrugMatches = textscan(fid, '%s%s','Delimiter','\t');
fclose(fid);

gene_drug_mapping = struct;
gene_drug_mapping.genes = geneDrugMatches{1};
gene_drug_mapping.drugIDs = geneDrugMatches{2};

% Export
StrepGeneDrugPairs = struct;
StrepGeneDrugPairs.mitis_gene_drug_mapping = gene_drug_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus gallolyticus strain ICDDRB-NRC-S3 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('s.gallolyticus.results.targets','r');
geneDrugMatches = textscan(fid, '%s%s','Delimiter','\t');
fclose(fid);

gene_drug_mapping = struct;
gene_drug_mapping.genes = geneDrugMatches{1};
gene_drug_mapping.drugIDs = geneDrugMatches{2};

% Export
StrepGeneDrugPairs.gallolyticus_gene_drug_mapping = gene_drug_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus oralis ATCC 49296 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('s.oralis.results.targets','r');
geneDrugMatches = textscan(fid, '%s%s','Delimiter','\t');
fclose(fid);

gene_drug_mapping = struct;
gene_drug_mapping.genes = geneDrugMatches{1};
gene_drug_mapping.drugIDs = geneDrugMatches{2};

% Export
StrepGeneDrugPairs.oralis_gene_drug_mapping = gene_drug_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus equinus strain AG46 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('s.equinus.results.targets','r');
geneDrugMatches = textscan(fid, '%s%s','Delimiter','\t');
fclose(fid);

gene_drug_mapping = struct;
gene_drug_mapping.genes = geneDrugMatches{1};
gene_drug_mapping.drugIDs = geneDrugMatches{2};

% Export
StrepGeneDrugPairs.equinus_gene_drug_mapping = gene_drug_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus pneumoniae genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('s.pneumoniae.results.targets','r');
geneDrugMatches = textscan(fid, '%s%s','Delimiter','\t');
fclose(fid);

gene_drug_mapping = struct;
gene_drug_mapping.genes = geneDrugMatches{1};
gene_drug_mapping.drugIDs = geneDrugMatches{2};

% Export
StrepGeneDrugPairs.pneumoniae_gene_drug_mapping = gene_drug_mapping;

%------------------------------------------------------------------------
% Import reactions from Streptococcus vestibularis strain 22-06 S6 genome
%------------------------------------------------------------------------
% Schema: ID	Name	Equation	Definition	Genes
fid = fopen('s.vestibularis.results.targets','r');
geneDrugMatches = textscan(fid, '%s%s','Delimiter','\t');
fclose(fid);

gene_drug_mapping = struct;
gene_drug_mapping.genes = geneDrugMatches{1};
gene_drug_mapping.drugIDs = geneDrugMatches{2};

% Export
StrepGeneDrugPairs.vestibularis_gene_drug_mapping = gene_drug_mapping;

end