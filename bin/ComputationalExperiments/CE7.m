% Computational Experiment
% Predict gene essentiality 
%
% Written by Matt Biggs, 2016

ensembleList = {'CE7_ensemble_25gcs'};

n_gcs_list = [2,5,10,15,20,25,30];            
            
% Load ensembles
if ~exist('CE7_ensemble_25gcs','var')
    load CE7_ensemble_25gcs
    CE7_ensemble_25gcs(cellfun(@isempty,CE7_ensemble_25gcs)) = [];
end

% Evaluate the networks individually and as ensembles
N = 21;
networkAccuracy = zeros(N*7,2);
ensembleAccuracy = zeros(7,2);
ensemblePrecisionsNRecalls = zeros(0,4);
for i = 1:length(n_gcs_list)
    curEnsemble = eval(ensembleList{i});
    gcs = n_gcs_list(i);
    N = length(curEnsemble);
    
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
    
    % Evaluate ensemble accuracies
    testGCs = zeros(length(testGCs),N);
    testNGCs = zeros(length(testNGCs),N);
    for j = 1:N
        testGCs(:,j) = curEnsemble{j}.gc_bm_vals(31:end);
        testNGCs(:,j) = curEnsemble{j}.ngc_bm_vals(1:17);
    end
%     testGCs = sum(testGCs > 1e-10,2) >= 11;
%     testNGCs = sum(testNGCs > 1e-10,2) >= 11;
    testGCs = sum(testGCs > 1e-10,2) >= N/2;
    testNGCs = sum(testNGCs > 1e-10,2) >= N/2;
    TP = sum( testGCs );
    TN = sum( testNGCs == 0 );
    FP = sum( testNGCs );
    FN = sum( testGCs == 0 );
    ensembleAccuracy(i,1) = (TP + TN) / (TP + TN + FP + FN);
    ensembleAccuracy(i,2) = gcs;
    
    % Evaluate ensemble precision/recall
    for t = 1:N
        testGCs = zeros(length(testGCs),N);
        testNGCs = zeros(length(testNGCs),N);
        for j = 1:N
            testGCs(:,j) = curEnsemble{j}.gc_bm_vals(31:end);
            testNGCs(:,j) = curEnsemble{j}.ngc_bm_vals(1:17);
        end
        testGCs = sum(testGCs > 1e-10,2) >= t;
        testNGCs = sum(testNGCs > 1e-10,2) >= t;
        TP = sum( testGCs );
        TN = sum( testNGCs == 0 );
        FP = sum( testNGCs );
        FN = sum( testGCs == 0 );
        ensemblePrecisionsNRecalls(end+1,1) = TP / (TP + FP); % precision
        ensemblePrecisionsNRecalls(end,2) = TP / (TP + FN); % recall
        ensemblePrecisionsNRecalls(end,3) = t;
        ensemblePrecisionsNRecalls(end,4) = gcs;
    end
end

% Write to file
% dlmwrite('CE7_networkAccuracies.tsv',networkAccuracy,'\t');
% dlmwrite('CE7_ensembleAccuracy.tsv',ensembleAccuracy,'\t');
% dlmwrite('CE7_ensemblePrecisionsNRecalls.tsv',ensemblePrecisionsNRecalls,'\t');
% 
% 
