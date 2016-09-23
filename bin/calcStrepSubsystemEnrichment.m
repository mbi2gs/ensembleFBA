function [subsysPvals,uSubsys] = calcStrepSubsystemEnrichment(inFilePrefix,N)
%---------------------------------------------------------------------
% calcStrepSubsystemEnrichment - Count the essential and non-essential rxns
% and calculate enrichment scores based on the rxn subsystems.
%
% Written by Matt Biggs, 2016
%---------------------------------------------------------------------

uRxns = {};
for i = 1:N
    % Read in file
    fid = fopen([inFilePrefix num2str(i) '.tsv'],'r');
    rxns = textscan(fid, '%s%d','Delimiter','\t');
    fclose(fid);
    
    uRxns = [uRxns; rxns{1}];
end
uRxns = unique(uRxns);

% Get the rxn essentiality information
rxnEssentialityMat = zeros(length(uRxns),N);
for i = 1:N
    % Read in file
    fid = fopen([inFilePrefix num2str(i) '.tsv'],'r');
    rxns = textscan(fid, '%s%d','Delimiter','\t');
    fclose(fid);
    
    mappedEssentiality = zeros(size(uRxns));
    for j = 1:length(rxns{1})
        curRxn = rxns{1}{j};
        currEss = rxns{2}(j);        
        mappedEssentiality(ismember(uRxns,curRxn)) = currEss;
    end
    rxnEssentialityMat(:,i) = mappedEssentiality;
end

essentialRxnIndicators = sum(rxnEssentialityMat,2) > N/2;
essentialRxns = uRxns(essentialRxnIndicators);

% Read in subsystem-to-reaction mapping
fid = fopen('SEED.rxn.pathways.tsv','r');
rxnSubsysMap = textscan(fid, '%s%s','Delimiter','\t');
fclose(fid);

uSubsys = {};
for i = 1:length(rxnSubsysMap{2})
    subsysList = strsplit(rxnSubsysMap{2}{i},'|');
    for j = 1:length(subsysList)
        uSubsys = [uSubsys; subsysList(j)];
    end
end
uSubsys = unique(uSubsys);

% Calculate enrichment scores using hypergeometric distribution
M = length(uRxns);      % "population size"
N = sum(essentialRxnIndicators); % "number of samples drawn"

ks = zeros(size(uSubsys));
xs = zeros(size(uSubsys));
subsysPvals = zeros(size(uSubsys));
for i = 1:length(uSubsys)
    % How many total reactions are annotated with each subsys
    uRxnsSubsystems = rxnSubsysMap{2}(ismember(rxnSubsysMap{1},uRxns));
    ks(i) = sum(~cellfun(@isempty,strfind(uRxnsSubsystems,uSubsys{i})));
    
    % How many essential reactions are annotated with each subsys
    essRxnsSubsystems = rxnSubsysMap{2}(ismember(rxnSubsysMap{1},essentialRxns));
    xs(i) = sum(~cellfun(@isempty,strfind(essRxnsSubsystems,uSubsys{i})));
    
    % Calculate enrichment (smaller p-value = greater enrichment)
%     p = hygepdf(x:k,M,k,N)
%     subsysPvals(i) = p;
end

dlmwrite([inFilePrefix 'enrichment.tsv'], [ [M N] ; [ks xs] ], '\t');







end