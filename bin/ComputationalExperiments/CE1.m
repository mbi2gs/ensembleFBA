% Computational Experiment
% Show that order matters when gap filling sequentially
%
% Written by Matt Biggs, 2016

% Load universal reaction database and add exchange rxns
load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Get the PA14 data formatted with work with the SEED database
[PA14Data] = getPA14GrowthConditions(seed_rxns_mat);

% Get the PA14 data formatted with work with the SEED database
%     PA14Data.biomassFn,growthCarbonSources,growthConditions,nonGrowthCarbonSources,nonGrowthConditions
[PA14Data] = getPA14GrowthConditions(seed_rxns_mat);

% 
removeBlockedBiomassComponents(seed_rxns_mat,PA14Data.growthConditions(:,1), PA14Data.biomassFn)
