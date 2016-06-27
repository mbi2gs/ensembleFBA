function [gc_bm_vals,ngc_bm_vals] = testModelsInGrowthConditions_flex(modelList,exchangeRxnNameList,growthConditions,nonGrowthConditions,verbose)
% Check "growth" (flux through biomass obj using FBA) under each set of conditions
% 
% Written by Matt Biggs, UVA, 2016

if verbose > 0
    fprintf('Testing each model in all growth and nongrowth conditions\n');
end

gc_bm_vals = zeros(size(growthConditions,2),length(modelList));
ngc_bm_vals = zeros(size(nonGrowthConditions,2),length(modelList));

for i = 1:length(modelList)
    % Calculate
    for j = 1:size(growthConditions,2)
        gc_bm_vals(j,i) = fba_flex(modelList{i},exchangeRxnNameList,growthConditions(:,j),1);
    end

    for j = 1:size(nonGrowthConditions,2)
        ngc_bm_vals(j,i) = fba_flex(modelList{i},exchangeRxnNameList,nonGrowthConditions(:,j),verbose);
    end
end

if verbose > 0
    fprintf('[Conditions (rows), Models (cols)]\n');
    gc_bm_vals
    ngc_bm_vals
end

end