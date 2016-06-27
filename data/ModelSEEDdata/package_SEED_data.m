% Read in SEED reactions and make a COBRA-format object
% with all reactions and compounds

fid = fopen('compound_info_SEED.tsv','r');
cmpds = textscan(fid, '%s%s%s','Delimiter','\t');
fclose(fid);

fid = fopen('complete_SEED_matrix.tsv','r');
rxns = textscan(fid, '%s%d%s%s','Delimiter','\t');
fclose(fid);

% Creat COBRA-format object
seed_rxns_mat.mets = cmpds{1};
seed_rxns_mat.metFormulas = cmpds{2};
seed_rxns_mat.metNames = cmpds{3};
seed_rxns_mat.rxns = rxns{1};
seed_rxns_mat.rev = rxns{2};
seed_rxns_mat.rxnNames = rxns{4};
seed_rxns_mat.S = sparse(length(seed_rxns_mat.mets),length(seed_rxns_mat.rxns));

% Populate S matrix
for i = 1:length(rxns{3})
    line = rxns{3}{i};
    line = strrep(line, '(', '');
    line = strrep(line, ')', '');
    rps = strsplit(line,'|');
    for j = 1:length(rps)
       curRp = rps{j};
       if length(curRp) > 1
           entry = strsplit(curRp,',');
           compoundIndex = str2num(entry{1});
           coefficient = str2num(entry{2});
           seed_rxns_mat.S(compoundIndex, i) = coefficient;
       end
    end
end

save('seed_rxns.mat','seed_rxns_mat');

% Load in the EC-number-to-reaction mappings
fid = fopen('Enzyme_Class.aliases.txt','r');
ECaliases = textscan(fid, '%s%s%s','Delimiter','\t'); % EC \t Rxns \t Plant_stuff
fclose(fid);

ECrxns = {};
rxnECNumbers = {};
for i = 2:length(ECaliases{1}) % Start from 2 to skip the file header
    rxnList = strsplit(ECaliases{2}{i,:},'|');
    for j = 1:length(rxnList)
        ECrxns{end+1,1} = rxnList{j};
        rxnECNumbers{end+1,1} = ECaliases{1}{i,:};
    end
end

% Load in KEGG-ID-to-reaction mappings
fid = fopen('KEGG.aliases.txt','r');
KEGGaliases = textscan(fid, '%s%s%s','Delimiter','\t'); % KEGG \t Rxns \t Plant_stuff
fclose(fid);

KEGGrxns = {};
rxnKEGGids = {};
for i = 2:length(KEGGaliases{1}) % Start from 2 to skip the file header
    rxnList = strsplit(KEGGaliases{2}{i,:},'|');
    for j = 1:length(rxnList)
        KEGGrxns{end+1,1} = rxnList{j};
        rxnKEGGids{end+1,1} = KEGGaliases{1}{i,:};
    end
end

% Load in Rxn-Name-to-reaction mappings
fid = fopen('name.aliases.txt','r');
Namealiases = textscan(fid, '%s%s%s','Delimiter','\t'); % Name \t Rxns \t Plant_stuff
fclose(fid);

NAMErxns = {};
rxnNames = {};
for i = 2:length(Namealiases{1}) % Start from 2 to skip the file header
    rxnList = strsplit(Namealiases{2}{i,:},'|');
    for j = 1:length(rxnList)
        NAMErxns{end+1,1} = rxnList{j};
        rxnNames{end+1,1} = Namealiases{1}{i,:};
    end
end

% Create a list of EC numbers, Names and KEGG IDs which correspond to
% "seed_rxns_mat"
seed_aliases = struct;
seed_aliases.rxns = seed_rxns_mat.rxns;
seed_aliases.rxnNames = cell(size(seed_aliases.rxns));
seed_aliases.rxnECNumbers = cell(size(seed_aliases.rxns));
seed_aliases.rxnKEGGids = cell(size(seed_aliases.rxns));
for i = 1:length(seed_aliases.rxns)
    curRxn = seed_aliases.rxns{i};
    % Match Name
    matches = ismember(NAMErxns,curRxn);
    seed_aliases.rxnNames{i} = strjoin(rxnNames(matches)','|');
    
    % Match EC Num
    matches = ismember(ECrxns,curRxn);
    seed_aliases.rxnECNumbers{i} = strjoin(rxnECNumbers(matches)','|');
    
    % Match KEGG ID
    matches = ismember(KEGGrxns,curRxn);
    seed_aliases.rxnKEGGids{i} = strjoin(rxnKEGGids(matches)','|');
end

save('seed_aliases.mat','seed_aliases');

