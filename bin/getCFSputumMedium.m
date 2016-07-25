function [cfSputumMedium] = getCFSputumMedium(universalRxnSet)
%-------------------------------------------------------------------------- 
% getCFSputumMedium - Formats the growth conditions that mimic CF sputum in 
% a format that is compatible with the input rxn database.
%
% Inputs:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%       matrix for exchange rxns (X matrix), a reversability indicator for
%       all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange 
%       rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames), 
%       and metabolite formulas (metFormulas)
%
% Outputs:
%     cfSputumMedium - a column of lower bounds for the exchange rxns in universalRxnSet.X
%
% Written by Matt Biggs
%--------------------------------------------------------------------------

%## Define growth conditions (based on synthetic CF medium conditions from iPAU1129)
% H2O               cpd00001	7930   (Index)
% O2                cpd00007    7928
% Pi                cpd00009    7934
% CO2               cpd00011    2940
% NH4               cpd00013    2938
% L-Glutamate       cpd00023    7584
% Glucose           cpd00027    7580
% Glycine           cpd00033    2573
% D-Alanine         cpd00117    980
% L-Lysine          cpd00039    2567
% L-Aspartate       cpd00041    7223
% SO4               cpd00048    7222
% L-Arginine        cpd00051    3141
% L-Serine          cpd00054    2202
% Methionine        cpd00060    6861
% Ornithine         cpd00064    6865
% L-Tryptophan      cpd00065    6864
% L-Phenylalanine   cpd00066    6863
% H+                cpd00067    6862
% L-Tyrosine        cpd00069    6867
% L-Cysteine        cpd00084    6520
% L-Leucine         cpd00107    5622
% L-Histidine       cpd00119    978
% L-Proline         cpd00129    6008
% L-Valine          cpd00156    230
% L-Threonine       cpd00161    5266
%%% No reactions use: K+ cpd00205   
% Nitrate           cpd00209    2205
% L-Lactate         cpd00159    233
% Mg2               cpd00254    6753
% L-Isoleucine      cpd00322    4296
%%% No reactions use: Na        
% Fe3               cpd10516    3321

cfSputumMedium = zeros(length(universalRxnSet.Ex_names),1);
cfSputumMedium([7930,7928,7934,2940,2938,7584,7580,2573,980,2567, ...
				7223,7222,3141,2202,6861,6865,6864,6863,6862,6867, ...
				6520,5622,978,6008,230,5266,2205,233,6753,4296, ...
				3321],1) = -1000;

end