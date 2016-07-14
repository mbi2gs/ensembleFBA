% Computational Experiment
% Remove biomass components that are not producible given the growth
% conditions and the reactions in the SEED database
%
% Written by Matt Biggs, 2016

load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

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
% Valine                    cpd00156    230
% Isoleucine                cpd00322    4296
% Leucine                   cpd00107    5622
% Methionine                cpd00060    6861
% Phenylalanine             cpd00066    6863
% Tyrosine                  cpd00069    6867
% Tryptophan                cpd00065    6864
% Cardiolipin               cpd12801    6246
% CoA                       cpd00010    2939
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
% Phosphatidylethanolamine  cpd11456    6907
% peptidoglycan subunit     cpd16008    3356
% Phosphatidylglycerol      cpd11652    7958
% phosphatidylserine        cpd15657    6603
% Putrescine                cpd00118    979
% Spermidine                cpd00264    6545
% Succinyl-CoA              cpd00078    1865
% UDP-glucose               cpd00026    7581
% UTP                       cpd00062    6859
% 
% Products:
% ADP           cpd00008    7935
% H+            cpd00067    6862
% Pi            cpd00009    7934
% PPi           cpd00012    2937
biomassFn = zeros(length(seed_rxns_mat.mets),1);
biomassFn([7933,2573,978,2567,7223,7584,2202,5266,1349,3130, ...
            6520,6008,230,4296,5622,6861,6863,6867,6864,6246, ...
            2939,2204,982,635,7842,634,2944,231,2566,7930, ...
            7932,7927,7929,7926,6907,3356,7958,6603,979,6545, ...
            1865,7581,6859],1) = -1;
biomassFn([7935,6862,7934,2937],1) = 1;

%## Define growth conditions (based on minimal media conditions from iMO1056)
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
% 4-Hydroxybenzoate     cpd00136    1346
% Acetate               cpd00029    7578
% Citrate               cpd00137    1345
% Glycerol-3-phosphate  cpd00080    6524
% D-Alanine             cpd00117    980
% D-Fructose            cpd00082    6526
% Adenosine             cpd00182    4234
% Fumarate              cpd00106    5621
% L-Arginine            cpd00051    3141
% L-Asparagine          cpd00132    1349
% L-Aspartate           cpd00041    7223
% L-Glutamate           cpd00023    7584
% L-Glutamine           cpd00053    3130
% L-Histidine           cpd00119    978
% Glycerol              cpd00100    5623
% Glycine               cpd00033    2573
% Itaconate             cpd00380    8748 % Exclude
% Ornithine             cpd00064    6865
% L-Phenylalanine       cpd00066    6863
% L-Proline             cpd00129    6008
% L-Serine              cpd00054    2202
% Malonate              cpd00308    8195
% L-Alanine             cpd00035    2568
% Glucose               cpd00027    7580
% Putrescine            cpd00118    979
% Pyruvate              cpd00020    7586
% Succinate             cpd00036    2569
% 2-Oxoglutarate        cpd00024    7583
% GABA                  cpd00281    8288
% L-Isoleucine          cpd00322    4296
% L-Malate              cpd00130    1347
% L-Lactate             cpd00159    233
% L-Leucine             cpd00107    5622
% D-Mannitol            cpd00314    9684
% Aminoethanol          cpd00162    5263
% N-Acetyl-L-glutamate  cpd00477    3940
% gly-glu-L             cpd11592    3895
% gly-pro-L             cpd11588    8803
% Hydroxy-L-proline     cpd00851    1958
% Butyrate              cpd00211    1874
% L-alanylglycine       cpd11585    3450
% L-Lysine              cpd00039    2567
% Carnitine             cpd00266    6547
% GLCN                  cpd00222    7244
% Hydroxyphenylacetate  cpd00489    247
% Propionate            cpd00141    4918
% Uridine               cpd00249    6209
% Hydroxybutanoate      cpd00797    7608

minimalMediaBase = zeros(length(seed_rxns_mat.mets),1);
minimalMediaBase([4926,2213,3321,6862,7930,6753,6714,2938,7928,7934,7222],1) = -1000;

growthCarbonSources = [1346,7578,1345,6524,980,6526,4234,5621,3141, ...
                       1349,7223,7584,3130,978,5623,2573,6865,6863,6008, ...
                       2202,8195,2568,7580,979,7586,2569,7583,8288,4296,...
                       1347,233,5622,9684,5263,3940,3895,8803,1958,1874, ...
                       3450,2567,6547,7244,247,4918,6209,7608];
               
n = length(growthCarbonSources);
growthConditions = repmat(minimalMediaBase,[1,n]);
for i = 1:n
    growthConditions(growthCarbonSources(i),i) = -10;
end


[newBiomassFn] = removeBlockedBiomassComponents(seed_rxns_mat,growthConditions(:,2),biomassFn);

fbmf = find(biomassFn);
toRemove = fbmf(~ismember(fbmf,find(newBiomassFn)));

% Need to remove 7 components from the biomass function (can't be produced 
% on the first several carbon sources):
% Valine                    cpd00156    230 
% Succinyl-CoA              cpd00078    1865  
% CoA                       cpd00010    2939  
% Cardiolipin               cpd12801    6246  
% phosphatidylserine        cpd15657    6603  
% Phosphatidylethanolamine  cpd11456    6907  
% Phosphatidylglycerol      cpd11652    7958  

