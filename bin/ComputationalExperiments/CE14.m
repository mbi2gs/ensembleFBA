% Computational Experiment
% Evaluate ensembles where individual models are trained on varying
% percentages of the total data set
%
% Written by Matt Biggs, 2016

ensembleList = {'ensemble_0p2', ...
                'ensemble_0p4', ...
                'ensemble_0p6', ...
                'ensemble_0p8', ...
                'ensemble_1'};

% Load ensembles
for i = 1:length(ensembleList)
    if ~exist(ensembleList{i},'var')
        eval(['load CE14_' ensembleList{i}]);
        eval([ensembleList{i} '(cellfun(@isempty,' ensembleList{i} ')) = []']);
    end
end

% Evaluate the networks individually and as ensembles
fraction_list = [0.2,0.4,0.6,0.8,1];
N = 21;

ensemblePrecision = zeros(length(fraction_list),4);
ensembleRecall = zeros(length(fraction_list),4);
ensembleAccuracy = zeros(length(fraction_list),4);
for i = 1:length(fraction_list)
    curEnsemble = eval(ensembleList{i});
    frac = fraction_list(i);
 
    % Evaluate ensembles
    testGCs = zeros(length(curEnsemble{1}.gc_bm_vals(21:end)),N);
    testNGCs = zeros(length(curEnsemble{1}.ngc_bm_vals(11:end)),N);
    for j = 1:N
        testGCs(:,j) = curEnsemble{j}.gc_bm_vals(21:end);
        testNGCs(:,j) = curEnsemble{j}.ngc_bm_vals(11:end);
    end
    
    thresholds = [1,11,21];
    for t = 1:3
        thresh = thresholds(t);
        testGCs2 = sum(testGCs > 1e-10,2) >= thresh;
        testNGCs2 = sum(testNGCs > 1e-10,2) >= thresh;
        TP = sum( testGCs2 );
        TN = sum( testNGCs2 == 0 );
        FP = sum( testNGCs2 );
        FN = sum( testGCs2 == 0 );
        ensemblePrecision(i,t) = TP / (TP + FP);
        ensembleRecall(i,t) = TP / (TP + FN);
        ensembleAccuracy(i,t) = (TP + TN) / (TP + TN + FP + FN);
    end
    ensemblePrecision(i,4) = frac;
    ensembleRecall(i,4) = frac;
    ensembleAccuracy(i,4) = frac;
end

% Write to file
dlmwrite('CE14_ensemblePrecision.tsv',ensemblePrecision,'\t');
dlmwrite('CE14_ensembleRecall.tsv',ensembleRecall,'\t');
dlmwrite('CE14_ensembleAccuracy.tsv',ensembleAccuracy,'\t');

