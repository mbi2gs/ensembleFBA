function [StrepData] = getStrepGrowthConditions(universalRxnSet)
%-------------------------------------------------------------------------- 
% getStrepGrowthConditions - Generates a biomass function, and formats the
% growth conditions for use with the input rxn database.
%
% Inputs:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%       matrix for exchange rxns (X matrix), a reversability indicator for
%       all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange 
%       rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames), 
%       and metabolite formulas (metFormulas)
%
% Outputs:
%     StrepData - a Matlab struct with the following fields:
%       biomassFn - the same format as a rxn (column) in universalRxnSet.S
%       growthCarbonSources - the list of carbon sources which allow growth 
%       growthConditions - a matrix of lower bounds for the exchange rxns in universalRxnSet.X
%       growthIndicators - a binary matrix indicating the carbon sources in
%                          which each species can grow (1 = growth)
%       speciesOrder - a cell array of species names corresponding to
%                      columns in "growthIndicators".
%
% Written by Matt Biggs
%--------------------------------------------------------------------------

%## Build biomass function (based on iMO1056):
% Substrates:
% ATP                       cpd00002	7933 (index)
% Glycine                   cpd00033	2573
% Histidine                 cpd00119	978
% Lysine                    cpd00039    2567
% Aspartate                 cpd00041    7223
% Glutamate                 cpd00023    7584
% Serine                    cpd00054    2202
% Threonine                 cpd00161    5266
% Asparagine                cpd00132    1349
% Glutamine                 cpd00053    3130
% Cysteine                  cpd00084    6520
% Proline                   cpd00129    6008
% Isoleucine                cpd00322    4296
% Leucine                   cpd00107    5622
% Methionine                cpd00060    6861
% Phenylalanine             cpd00066    6863
% Tyrosine                  cpd00069    6867
% Tryptophan                cpd00065    6864
% CTP                       cpd00052    2204
% dATP                      cpd00115    982
% dCTP                      cpd00356    635
% dGTP                      cpd00241    7842
% dTTP                      cpd00357    634
% FAD                       cpd00015    2944
% Glycogen                  cpd00155    231
% GTP                       cpd00038    2566
% H2O                       cpd00001    7930
% NAD                       cpd00003    7932
% NADH                      cpd00004    7927
% NADP                      cpd00006    7929
% NADPH                     cpd00005    7926
% peptidoglycan subunit     cpd16008    3356
% Putrescine                cpd00118    979
% Spermidine                cpd00264    6545
% UDP-glucose               cpd00026    7581
% UTP                       cpd00062    6859
% 
% Products:
% ADP           cpd00008    7935
% H+            cpd00067    6862
% Pi            cpd00009    7934
% PPi           cpd00012    2937
biomassFn = zeros(length(universalRxnSet.Ex_names),1);
biomassFn([7933,2573,978,2567,7223,7584,2202,5266,1349,3130, ...
            6520,6008,4296,5622,6861,6863,6867,6864, ...
            2204,982,635,7842,634,2944,231,2566,7930, ...
            7932,7927,7929,7926,3356,979,6545, ...
            7581,6859],1) = -1;
biomassFn([7935,6862,7934,2937],1) = 1;

%## Define growth conditions
% Cobalt        cpd00149    4926 (Index)
% Copper        cpd00058    2213
% Fe3           cpd10516    3321
% H+            cpd00067    6862
% H2O           cpd00001    7930
% Mg2           cpd00254    6753
% N2            cpd00528    6714
% NH4           cpd00013    2938
% O2            cpd00007    7928
% Pi            cpd00009    7934
% SO4           cpd00048    7222
%%% No reactions use: Cd2+          cpd01012     
%%% No reactions use: Fe2           cpd00021    
%%% No reactions use: K+            cpd00205    
%%% No reactions use: Mn2           cpd00030    
%%% No reactions use: Na            cpd00971    
%%% No reactions use: Zn2           cpd00034    
% Growth Carbon sources:
% Dextrin                       cpd11594    8456
% Maltose                       cpd00179    602
% D-cellobiose                  cpd03844    5204
% Gentiobiose                   cpd05158    9463
% Sucrose                       cpd00076    1863
% Stachyose                     cpd01133    4431
% a-D-Lactose                   cpd01354    1766
% D-Melibiose                   cpd03198    1783
% Salicin                       cpd01030    1884
% N-Acetyl-DGlucosamine         cpd00122    5998
% N-Acetyl-b-D-Mannosamine      cpd00492    4947
% N-Acetyl-Neuraminic Acid      cpd27569    3485
% a-D-Glucose                   cpd00027    7580
% D-Mannose                     cpd00138    1351
% D-Fructose                    cpd00082    6526
% D-Galactose                   cpd00709    1202
% L-Fucose                      cpd00751    5461
% L-Rhamnose                    cpd00396    1368
% Inosine                       cpd00246    6219
% D-Sorbitol                    cpd00588    7827
% D-Mannitol                    cpd00314    9684
% M-Inositol                    cpd00121    6001
% Glycerol                      cpd00100    5623
% Glucose-6-Phosphate           cpd00079    1866
% Fructose-6-Phosphate          cpd00072    1859
% D-Aspartic Acid               cpd00320    4294
% D-Serine                      cpd00550    2431
% L-Alanine                     cpd00035    2568
% L-Arginine                    cpd00051    3141
% L-Aspartic Acid               cpd00041    7223
% L-Glutamic Acid               cpd00023    7584
% L-Histidine                   cpd00119    978
% L-Serine                      cpd00054    2202
% D-GalacturonicAcid            cpd00280    8289
% D-GlucuronicAcid              cpd00164    5261
% Quinic Acid                   cpd00248    6210
% Acetic Acid                    cpd00029    7578
% L-Lactic Acid                 cpd00159    233
% Citric Acid                   cpd00137    1345
% D-Malic Acid                  cpd00386    6024
% L-Malic Acid                  cpd00130    1347
% Acetoacetic Acid              cpd00142    4919
% Propionic Acid                cpd00141    4918

minimalMediaBase = zeros(length(universalRxnSet.Ex_names),1);
minimalMediaBase([4926,2213,3321,6862,7930,6753,6714,2938,7928,7934,7222],1) = -1000;

growthCarbonSources = [8456,602,5204,9463,1863,4431,1766,1783,1884,5998, ...
                       4947,3485,7580,1351,6526,1202,5461,1368,6219,7827, ...
					   9684,6001,5623,1866,1859,4294,2431,2568,3141,7223, ...
					   7584,978,2202,8289,5261,6210,7578,233,1345, ...
					   6024,1347,4919,4918 ];

n = length(growthCarbonSources);
growthConditions = repmat(minimalMediaBase,[1,n]);
for i = 1:n
    growthConditions(growthCarbonSources(i),i) = -10;
end

StrepData = struct;
StrepData.biomassFn = biomassFn;
StrepData.minimalMediaBase = minimalMediaBase;
StrepData.growthCarbonSources = growthCarbonSources;
StrepData.growthConditions = growthConditions;

% Indicate which species grow on which carbon sources
growthIndicators = [52	78	6	0	80	49	70	28	0	72	70	61	79	60	73	65	0	0	2	0	0	0	19	13	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	25	4	19; ...
                    81	78	78	79	75	77	77	26	77	79	0	0	77	79	77	79	0	0	0	0	78	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	31	3	12; ...
                    24	59	5	0	62	0	34	0	0	49	15	0	54	53	56	35	0	0	0	1	0	0	12	2	0	0	0	0	0	0	0	0	6	0	0	9	0	0	0	0	49	12	37; ...
                    57	74	66	67	74	55	70	63	58	76	17	3	52	66	72	59	0	0	0	0	0	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	37	4	12; ...
                    60	78	16	35	72	30	76	19	23	76	76	58	74	71	70	71	0	0	4	0	0	0	40	0	0	0	0	0	0	0	0	0	1	0	0	13	0	0	0	0	53	7	19; ...
                    67  79	82	82	81	0	78	0	82	27	0	0	82	80	80	78	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	4	0	0	0	0	0	0	0	0	24	5	11];

growthIndicators(growthIndicators > 0) = 1;                
StrepData.growthIndicators = growthIndicators';
StrepData.speciesOrder = {'S. mitis', 'S. gallolyticus', 'S. oralis', 'S. equinus', 'S. pneumoniae', 'S. vestibularis'};

end