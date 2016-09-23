% Computational Experiment
% Reconstruct and analyze ensembles for 6 Streptococcus species
%
% Written by Matt Biggs, 2016

% Load universal reaction database and add exchange rxns
load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Get the Streptococcus data formatted to work with the SEED database
[StrepData] = getStrepGrowthConditions(seed_rxns_mat);

% Get the Streptococcus gene-to-reaction mappings
[StrepGenomicData] = getStrepGenomeAnnotations();
rxnMappingsList = fieldnames(StrepGenomicData);

% Get the Streptococcus gene/drug mappings
[StrepGeneDrugPairs] = getStrepGeneDrugPairs();

N = 21;
ensembles = cell(N,6);
for i = 1:length(rxnMappingsList)
    fprintf(['\titeration ' num2str(i) '\n']);
    rxnMappingsList{i}
    
    % Force networks to contain reactions annotated from the genome
    rxn_GPR_mapping = getfield(StrepGenomicData,rxnMappingsList{i});
    Urxns2set = [find(ismember(seed_rxns_mat.rxns,rxn_GPR_mapping.rxns)); find(ismember(seed_rxns_mat.rxns,'rxn05064'))]; % include spontaneous rxn05064
    Uset2 = ones(size(Urxns2set));
    
    % Include exchange reactions for all non-growth conditions, just so that
    % it's the network itself--not the lack of exchange reactions--that prevents growth
    Xrxns2set = find(sum( abs(seed_rxns_mat.X(StrepData.growthCarbonSources(:),:)) ,1) > 0);
    Xset2 = ones(size(Xrxns2set));
    
    % Set parameters
    biologicalData = struct;
    biologicalData.biomassFn = StrepData.biomassFn;
    biologicalData.Urxns2set = Urxns2set;
    biologicalData.Uset2 = Uset2;
    biologicalData.Xrxns2set = Xrxns2set;
    biologicalData.Xset2 = Xset2;
    
    % Randomly select growth conditions
    currGrowthIndicators = StrepData.growthIndicators(:,i);
    biologicalData.growthConditions = StrepData.growthConditions(:,currGrowthIndicators == 1);
    biologicalData.nonGrowthConditions = [];
    
    params = struct;
    params.sequential = 1;
    params.stochast = 1;
    params.numModels2gen = N;
    params.fractionUrxns2set = 0.8;
    params.rndSeed = i;
    params.numGrowthConditions = 25;
    params.numNonGrowthConditions = 0;
    params.verbose = 0;
    
    tic
    [modelList] = build_ensemble(seed_rxns_mat,biologicalData,params);
    stseq1 = toc;
    
    modelList2 = addGPRs(modelList,rxn_GPR_mapping);
    ensembles(:,i) = modelList2;
end
save('CE13_streptococcus_ensembles.mat','ensembles');

% Write ensemble to separate files for upload to HPC
% for i = 1:6
%    for j = 1:21 
%        m = ensembles{j,i};
%        fileName = [StrepData.speciesOrder{i} '_' num2str(j) '.mat']; 
%        fileName = strrep(fileName,'S. ','S.');
%        save(fileName,'m');
%    end
% end

load CE13_streptococcus_ensembles.mat

% Evaluate essential genes for all species and find overlap in drugs which
% will influence them
richMedia = sum(StrepData.growthConditions,2);
richMedia(richMedia < -1000) = -1000;
geneMappingsList = fieldnames(StrepGeneDrugPairs);
for i = 1:length(rxnMappingsList)
    curEnsemble = ensembles(:,i);
    
    % Identify genes that correspond to drug targets that are also
    % represented in the models
    allGenes = cell(0,1);
    for j = 1:N
        allGenes = [allGenes; curEnsemble{j}.genes];
    end
    allGenes = unique(allGenes);
    
    genesWdrugMatches = getfield(StrepGeneDrugPairs,geneMappingsList{i});
    genesInModelWdrugMatches = allGenes(ismember(allGenes,genesWdrugMatches.genes));
    
    % Check gene essentiality
    geneEssentialityByNet = zeros(length(genesInModelWdrugMatches),N);
    for j = 1:N
        curMod = curEnsemble{j};
        for k = 1:length(genesInModelWdrugMatches)
            curGene = genesInModelWdrugMatches{k};
            delMod = simulateGeneDeletion(curMod,curGene);
            delGrowth = fba_flex(delMod,seed_rxns_mat.Ex_names,richMedia,0);
            geneEssentialityByNet(k,j) = delGrowth < 1e-10;
        end
    end
    
    % Decide on essential genes
    essentialGeneIndicators = sum(geneEssentialityByNet,2) > N/2;
    essentialGenes = genesInModelWdrugMatches(essentialGeneIndicators > 0);
    drugsHitEssentialGenes = genesWdrugMatches.drugIDs(ismember(genesWdrugMatches.genes,essentialGenes));
    
    fileName = ['CE13_geneEssentiality_' StrepData.speciesOrder{i} '.mat']; fileName = strrep(fileName,'S. ','');
    save(fileName,'geneEssentialityByNet','genesInModelWdrugMatches','drugsHitEssentialGenes');
end

% Evaluate the drugs that uniquely interact with each species (and the
% genes they target)
load CE13_geneEssentiality_mitis.mat
mitis_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       mitis_drugs = [mitis_drugs; curDrug];
   end
end
mitis_drugs = unique(mitis_drugs);

load CE13_geneEssentiality_gallolyticus.mat
gallolyticus_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       gallolyticus_drugs = [gallolyticus_drugs; curDrug];
   end
end
gallolyticus_drugs = unique(gallolyticus_drugs);

load CE13_geneEssentiality_oralis.mat
oralis_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       oralis_drugs = [oralis_drugs; curDrug];
   end
end
oralis_drugs = unique(oralis_drugs);

load CE13_geneEssentiality_equinus.mat
equinus_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       equinus_drugs = [equinus_drugs; curDrug];
   end
end
equinus_drugs = unique(equinus_drugs);

load CE13_geneEssentiality_pneumoniae.mat
pneumoniae_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       pneumoniae_drugs = [pneumoniae_drugs; curDrug];
   end
end
pneumoniae_drugs = unique(pneumoniae_drugs);

load CE13_geneEssentiality_vestibularis.mat
vestibularis_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       vestibularis_drugs = [vestibularis_drugs; curDrug];
   end
end
vestibularis_drugs = unique(vestibularis_drugs);

size(mitis_drugs)
size(gallolyticus_drugs)
size(oralis_drugs)
size(equinus_drugs)
size(pneumoniae_drugs)
size(vestibularis_drugs)

unique2mitis = mitis_drugs(~ismember(mitis_drugs,[gallolyticus_drugs;oralis_drugs;equinus_drugs;pneumoniae_drugs;vestibularis_drugs]))
unique2gallolyticus = gallolyticus_drugs(~ismember(gallolyticus_drugs,[mitis_drugs;oralis_drugs;equinus_drugs;pneumoniae_drugs;vestibularis_drugs]))
unique2oralis = oralis_drugs(~ismember(oralis_drugs,[mitis_drugs;gallolyticus_drugs;equinus_drugs;pneumoniae_drugs;vestibularis_drugs]))
unique2equinus = equinus_drugs(~ismember(equinus_drugs,[mitis_drugs;gallolyticus_drugs;oralis_drugs;pneumoniae_drugs;vestibularis_drugs]))
unique2pneumoniae = pneumoniae_drugs(~ismember(pneumoniae_drugs,[mitis_drugs;gallolyticus_drugs;oralis_drugs;equinus_drugs;vestibularis_drugs]))
unique2vestibularis = vestibularis_drugs(~ismember(vestibularis_drugs,[mitis_drugs;gallolyticus_drugs;oralis_drugs;equinus_drugs;pneumoniae_drugs]))

allDrugs = [mitis_drugs;gallolyticus_drugs;oralis_drugs;equinus_drugs;pneumoniae_drugs;vestibularis_drugs];
uDrugs = unique(allDrugs);
counts = zeros(size(uDrugs));
for i = 1:length(uDrugs)
    counts(i) = sum(ismember(allDrugs,uDrugs{i}));
end
common2all = uDrugs(counts == 6)

%---------------------------------------------------------------------
% Evaluate the essential reactions in each species and the subsystem
% enrichment in each
%---------------------------------------------------------------------
% S. mitis
inputFile = 'EssRxns_S.mitis_';
[mitis_subsysPvals,uSubsys] = calcStrepSubsystemEnrichment(inputFile,N);

subsysEnrichment = zeros(length(uSubsys),6);
subsysEnrichment(:,1) = mitis_subsysPvals;

fid = fopen('CE13_uniqueSubsystems.txt','w');
for i = 1:length(uSubsys)
   fprintf(fid,[uSubsys{i} '\n']); 
end
fclose(fid);

% S. gallolyticus
inputFile = 'EssRxns_S.gallolyticus_';
[gallolyticus_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,2) = gallolyticus_subsysPvals;

% S. oralis
inputFile = 'EssRxns_S.oralis_';
[oralis_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,3) = oralis_subsysPvals;

% S. equinus
inputFile = 'EssRxns_S.equinus_';
[equinus_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,4) = equinus_subsysPvals;

% S. pneumoniae
inputFile = 'EssRxns_S.pneumoniae_';
[pneumoniae_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,5) = pneumoniae_subsysPvals;

% S. vestibularis
inputFile = 'EssRxns_S.vestibularis_';
[vestibularis_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,6) = vestibularis_subsysPvals;


