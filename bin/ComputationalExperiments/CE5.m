% Computational Experiment
% Evaluate ensemble method against individual models
%
% Written by Matt Biggs, 2016

ensembleList = {'CE5_ensemble_2gcs', ...
                'CE5_ensemble_5gcs', ...
                'CE5_ensemble_10gcs', ...
                'CE5_ensemble_15gcs', ...
                'CE5_ensemble_20gcs', ...
                'CE5_ensemble_25gcs', ...
                'CE5_ensemble_30gcs'};

n_gcs_list = [2,5,10,15,20,25,30];            
            
% Load ensembles
for i = 1:length(ensembleList)
    if ~exist(ensembleList{i},'var')
        eval(['load ' ensembleList{i}]);
    end
end

% Evaluate the networks individually and as ensembles
N = 21;
networkAccuracy = zeros(N*7,2);
ensembleAccuracy = zeros(7,2);
for i = 1:length(n_gcs_list)
    curEnsemble = eval(ensembleList{i});
    gcs = n_gcs_list(i);
    
    % Evaluate individual networks
    for j = 1:N
        testGCs = curEnsemble{j}.gc_bm_vals(31:end);
        testNGCs = curEnsemble{j}.ngc_bm_vals(1:17);
        TP = sum( testGCs > 1e-10 );
        TN = sum( testNGCs < 1e-10 );
        FP = sum( testNGCs > 1e-10 );
        FN = sum( testGCs < 1e-10 );
        networkAccuracy((i-1)*N+j,1) = (TP + TN) / (TP + TN + FP + FN);
        networkAccuracy((i-1)*N+j,2) = gcs;
    end
    
    % Evaluate ensembles
    testGCs = zeros(length(testGCs),N);
    testNGCs = zeros(length(testNGCs),N);
    for j = 1:N
        testGCs(:,j) = curEnsemble{j}.gc_bm_vals(31:end);
        testNGCs(:,j) = curEnsemble{j}.ngc_bm_vals(1:17);
    end
    testGCs = sum(testGCs > 1e-10,2) >= 11;
    testNGCs = sum(testNGCs > 1e-10,2) >= 11;
    TP = sum( testGCs );
    TN = sum( testNGCs == 0 );
    FP = sum( testNGCs );
    FN = sum( testGCs == 0 );
    ensembleAccuracy(i,1) = (TP + TN) / (TP + TN + FP + FN);
    ensembleAccuracy(i,2) = gcs;
end

% Write to file
dlmwrite('CE5_networkAccuracies.tsv',networkAccuracy,'\t');
dlmwrite('CE5_ensembleAccuracy.tsv',ensembleAccuracy,'\t');



