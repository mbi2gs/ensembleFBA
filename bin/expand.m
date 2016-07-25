function [model, rxnDatabase, fromUthresh, fromXthresh, feasible] = expand(universalRxnSet,growthConditions,nonGrowthConditions,biomassFn,Urxns2set,Uset2,Xrxns2set,Xset2,verbose,stochast,rndSeed)
%-------------------------------------------------------------------------- 
% expand - Gap fill algorithm that identifies a minimum set of reactions to
% add (from a universal rxn database such as SEED) in order to allow flux
% through the biomass function given a set of growth conditions. Uses a
% linear variable to constrain the absolute value of fluxes.
%
% Inputs:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%       matrix for exchange rxns (X matrix), a reversability indicator for
%       all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange 
%       rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames), 
%       and metabolite formulas (metFormulas).
%     growthConditions - set of lower bounds corresponding to growth media conditions
%     nonGrowthConditions - set of lower bounds corresponding to non-growth media conditions
%     biomassFn - same format as a reaction in universalRxnSet.S
%     Urxns2set - A list of rxn indices (in universalRxnSet.S) that are forced to be included or excluded
%     Uset2 - List of 1's or 0's (inclusion or exclusion)
%     Xrxns2set - Same as Urxns2set, but for exchange rxns (universalRxnSet.X)
%     Xset2 - Same as Uset2 but for exchange rxns (universalRxnSet.X)
%     verbose - Indicates verbose output (1 is some detail, 2 is all details)
%     stochast - 1 indicates stochastic weights on the reactions in universalRxnSet.S, 0 indicates deterministic weights
%     rndSeed - Allows the user to manually set the random seed
%
% Outputs:
%     model - A metabolic network reconstruction in COBRA format as a matlab struct
%     rxnDatabase - The model in the same format as universalRxnSet
%     fromUthresh - A logical vector indicating the reactions from universalRxnSet.S which were included in model
%     fromXthresh - A logical vector indicating the reactions from universalRxnSet.X which were included in model
%     feasible - 1 indicates successful model reconstruction. 0 indicates a solver error (perhaps an infeasible optimization problem).
%
% Written by Matt Biggs, mb3ad@virginia.edu, 2016
%-------------------------------------------------------------------------- 
if verbose > 0
    fprintf('\nExpansion Step\n');
end

if stochast > 0
    rng(rndSeed,'twister');
end

% Initialize
prevSolutionSize = size(universalRxnSet.S,2) + size(universalRxnSet.X,2);
curRxnSet = universalRxnSet;
curGrowthConditions = growthConditions;
curNonGrowthConditions = nonGrowthConditions;
curUrxns2set = Urxns2set;
curUset2 = Uset2;
curXrxns2set = Xrxns2set;
curXset2 = Xset2;
feasible = 1;
iterate = 1;
numIterationsCompleted = 0;

% Iterate to refine solution
while iterate > 0
   
    if verbose > 1
        fprintf('Expand step; trim a little more\n');
    end
        
    % Initialize outputs
    model = cell(0,1);
    rxnDatabase = cell(0,1);
    fromUthresh = [];
    fromXthresh = [];

    % Append biomass reaction to U
    Ubm = [curRxnSet.S biomassFn(:)];

    % Set X (exchange reaction) matrix
    X = curRxnSet.X;
    Ex_names = curRxnSet.Ex_names;

    % Identify total number of conditions
    n_gc = size(curGrowthConditions,2);
    n_ngc = size(curNonGrowthConditions,2);
    n_cond = n_gc;

    % Initialize variables
    [n_mets,n_Urxns] = size(Ubm);
    [~,n_Xrxns] = size(X);
    c_bm = zeros(n_Urxns+n_Xrxns,1); 
    c_bm(n_Urxns) = 1;
    n_UX = n_Urxns+n_Xrxns;

    % Set proper bounds
    ub_U = 1000*ones(n_Urxns,1);
    lb_U = [-1000*curRxnSet.rev; 0];
    ub_X = 1000*ones(n_Xrxns,1);
    
    % Build objective
    if stochast > 0
        objective = [zeros(1,n_UX*n_cond ) 1+rand(1,n_Urxns) ones(1,n_Xrxns)];
    else
        objective = [zeros(1,n_UX*n_cond ) ones(1,n_Urxns) ones(1,n_Xrxns)];
    end
    vtype = repmat('C',n_UX*n_cond + n_UX,1);

    %## Build matrices
    % Build UX diagonal matrix
    UX = sparse([Ubm X]);
    UX_list = repmat({UX},n_cond,1);
    UX_diag = blkdiag(UX_list{:}); % Repeates [Ubm X] along the diagonal n_cond times

    if verbose > 1
        fprintf('Built UX\n');
    end

    % Build upper and lower bound matrices
    lb_mat = sparse(n_UX*n_cond,n_UX*n_cond+n_UX);
    lb_rhs = zeros(n_UX*n_cond,1);
    for i = 1:n_cond
        lb_mat(((i-1)*n_UX+1):i*n_UX,:) = [repmat(sparse(n_UX,n_UX),[1,i-1]) speye(n_UX) repmat(sparse(n_UX,n_UX),[1,n_cond-i]) sparse(n_UX,n_UX)];
        lb_rhs(((i-1)*n_UX+1):i*n_UX,1) = [lb_U; curGrowthConditions(:,i)];
    end

    ub_mat = sparse(n_UX*n_cond,n_UX*n_cond+n_UX);
    for i = 1:n_cond
        ub_mat(((i-1)*n_UX+1):i*n_UX,:) = [repmat(sparse(n_UX,n_UX),[1,i-1]) speye(n_UX) repmat(sparse(n_UX,n_UX),[1,n_cond-i]) sparse(n_UX,n_UX)];
    end

    lb_mat_z = sparse(n_UX*n_cond,n_UX*n_cond+n_UX);
    for i = 1:n_cond
        lb_mat_z(((i-1)*n_UX+1):i*n_UX,:) = [repmat(sparse(n_UX,n_UX),[1,i-1]) speye(n_UX) repmat(sparse(n_UX,n_UX),[1,n_cond-i]) speye(n_UX,n_UX)];
    end

    ub_mat_z = sparse(n_UX*n_cond,n_UX*n_cond+n_UX);
    for i = 1:n_cond
        ub_mat_z(((i-1)*n_UX+1):i*n_UX,:) = [repmat(sparse(n_UX,n_UX),[1,i-1]) speye(n_UX) repmat(sparse(n_UX,n_UX),[1,n_cond-i]) -speye(n_UX,n_UX)];
    end

    lb_z = [sparse(n_UX,n_UX*n_cond) speye(n_UX)];

    if verbose > 1
        fprintf('Built LB/UB mats\n');
    end

    % Build biomass threshold matrices
    c_gc_list = repmat({c_bm'},n_gc,1);
    c_gc_diag = sparse(blkdiag(c_gc_list{:}));

    % Force specific z variables to always be on
    on_rxns = zeros(n_UX,1);
    on_rxns(curUrxns2set(curUset2 > 0)) = 1;
    on_rxns(curXrxns2set(curXset2 > 0)+n_Urxns) = 1;
    on_rxns_diag = spdiags(on_rxns,0,n_UX,n_UX);
    on_rxns_diag(sum(on_rxns_diag,2) == 0,:) = [];
    [orr,~] = size(on_rxns_diag);
    on_rxns_diag = [sparse(orr,n_UX*n_cond) on_rxns_diag];
    on_rxs_RHS = ones(orr,1);
    
    % Force specific z variables to always be off
    off_rxns = zeros(n_UX,1);
    off_rxns(curUrxns2set(curUset2 == 0)) = 1;
    off_rxns(curXrxns2set(curXset2 == 0)+n_Urxns) = 1;
    off_rxns_diag = spdiags(off_rxns,0,n_UX,n_UX);
    off_rxns_diag(sum(off_rxns_diag,2) == 0,:) = [];
    [orrf,~] = size(off_rxns_diag);
    off_rxns_diag = [sparse(orrf,n_UX*n_cond) off_rxns_diag];
    off_rxs_RHS = zeros(orrf,1);

    % Put it all together
    A =  [UX_diag sparse(n_mets*n_cond,n_UX) ; ...  %(1) Steady state assumption     
          lb_mat; ...                               %(2) Flux bounds
          ub_mat; ...                               % 
          lb_mat_z; ...                             %(3) Flux bounded by z
          ub_mat_z; ...                             % 
          lb_z; ...                                 %(4) z >= 0
          c_gc_diag sparse(n_gc,n_UX);...           %(5) Lower bound on biomass flux through growth conditions
          on_rxns_diag; ...                         %(6) Constrain some z values to be > 0
          off_rxns_diag; ...                        %    Constrain some z values to be = 0
          ]; 

    clear UX_diag lb_mat ub_mat lb_mat_z ub_mat_z lb_z c_gc_diag c_ngc_diag c_gc_list c_ngc_list 
    if verbose > 1
        fprintf('Built A\n');
    end

    RHS = [ zeros(n_mets*n_cond,1); ...             %(1)  
            lb_rhs; ...                             %(2)
            1000*ones(n_UX*n_cond,1); ...           %
            zeros(n_UX*n_cond,1); ...               %(3)
            zeros(n_UX*n_cond,1); ...               %
            zeros(n_UX,1); ...                      %(4)
            0.05*ones(n_gc,1);...                   %(5)
            on_rxs_RHS; ...                         %(6)
            off_rxs_RHS; ...                        %
            ];                           

    sense = [repmat('=',n_mets*n_cond,1); ...       %(1)    
             repmat('>',n_UX*n_cond,1); ...         %(2)
             repmat('<',n_UX*n_cond,1); ...         %
             repmat('>',n_UX*n_cond,1); ...         %(3)
             repmat('<',n_UX*n_cond,1); ...         %
             repmat('>',n_UX,1); ...                %(4)
             repmat('>',n_gc,1); ...                %(5)
             repmat('>',orr,1); ...                 %(6)
             repmat('<',orrf,1); ...                %
             ];    

    if verbose > 1
        fprintf('Built entire problem\n');     
    end

    if verbose > 1
        fprintf(['\nObjective size:   ' num2str(size(objective)) '\n']);
        fprintf(['Vtype size:       ' num2str(size(vtype)) '\n']);
        fprintf(['A size:           ' num2str(size(A)) '\n']);
        fprintf(['RHS size:         ' num2str(size(RHS)) '\n']);
        fprintf(['Sense size:       ' num2str(size(sense)) '\n\n']);
    end
    
    try
        % Optimize
        clear gurobi_model;
        gurobi_model.A = A;
        gurobi_model.obj = objective;
        gurobi_model.lb = [-1e30*ones(1,n_UX*n_cond) zeros(1,n_UX)];
        gurobi_model.rhs = RHS;
        gurobi_model.sense = sense;
        gurobi_model.vtype = vtype;
        gurobi_model.modelsense = 'min';
        
        clear params;
        params.outputflag = 0;
        result = gurobi(gurobi_model, params);
        
        if verbose > 1
            disp(result)
        end
        
        n_front = n_UX*n_cond;
        r_u = result.x(n_front+1 : n_front+n_Urxns);
        rulog = r_u > 1e-10;
        rulog(end) = [];
        r_x = result.x(n_front+n_Urxns+1 : n_front+n_Urxns+n_Xrxns);
        rxlog = r_x > 1e-10;
        numIterationsCompleted = numIterationsCompleted + 1;
        
        if sum(rulog) == 0
            feasible = 0;
        end

        % Make COBRA model of current solution
        newModel = struct;
        newModel.mets = curRxnSet.mets;
        newModel.metFormulas = curRxnSet.metFormulas;
        newModel.metNames = curRxnSet.metNames;
        newModel.rxns = [curRxnSet.rxns(rulog); Ex_names(rxlog); {'Biomass'}];
        newModel.rxnNames = [curRxnSet.rxnNames(rulog); Ex_names(rxlog); {'Biomass'}];
        newModel.rev = double([curRxnSet.rev(rulog); zeros(sum(rxlog),1); 0]);
        newModel.S = double([curRxnSet.S(:,rulog) X(:,rxlog) biomassFn(:)]);
        newModel.c = double(zeros(sum(rulog) + 1 + sum(rxlog),1)); newModel.c(end,1) = 1;
        newModel.b = double(zeros(length(curRxnSet.mets),1));
        newModel.ub = double([ub_U([rulog; 1>10]); ub_X(rxlog,1); 1000]);
        newModel.lb = double([lb_U([rulog; 1>10]); curGrowthConditions(rxlog,1); 0]);
        model = newModel;

        % Format as a reaction database for the next iteration or the trim step
        newRDB = struct;
        newRDB.mets = curRxnSet.mets;
        newRDB.metFormulas = curRxnSet.metFormulas;
        newRDB.metNames = curRxnSet.metNames;
        newRDB.rxns = curRxnSet.rxns(rulog);
        newRDB.rxnNames = curRxnSet.rxnNames(rulog);
        newRDB.rev = curRxnSet.rev(rulog);
        newRDB.S = curRxnSet.S(:,rulog);
        newRDB.X = X(:,rxlog);
        newRDB.Ex_names = Ex_names(rxlog);
        newRDB.growthConditions = curGrowthConditions(rxlog,:);
        if size(curNonGrowthConditions,2) > 0
            newRDB.nonGrowthConditions = curNonGrowthConditions(rxlog,:);
        else
            newRDB.nonGrowthConditions = curNonGrowthConditions;
        end
        rxnDatabase = newRDB;
        fromUthresh = rulog;
        fromXthresh = rxlog;
        
        % Set up for next iteration
        curSol = size(newRDB.S,2) + size(newRDB.X,2);
        if curSol < prevSolutionSize
            iterate = 1;
            prevSolutionSize = curSol;
            
            % Fix indices for rxns that need to be kept
            tmpUrxns2set = [];
            tmpUset2 = [];
            for crs = 1:length(curUrxns2set)
                rxn = curRxnSet.rxns{curUrxns2set(crs)};
                rxni = find(ismember(newRDB.rxns,rxn));
                if ~isempty(rxni)
                    tmpUrxns2set = [tmpUrxns2set; rxni];
                    tmpUset2 = [tmpUset2; curUset2(crs)];
                end
            end
            
            tmpXrxns2set = [];
            tmpXset2 = [];
            for crs = 1:length(curXrxns2set)
                rxn = curRxnSet.Ex_names{curXrxns2set(crs)};
                rxni = find(ismember(newRDB.Ex_names,rxn));
                if ~isempty(rxni)
                    tmpXrxns2set = [tmpXrxns2set; rxni];
                    tmpXset2 = [tmpXset2; curXset2(crs)];
                end
            end
            
            curUrxns2set = tmpUrxns2set;
            curUset2 = tmpUset2;
            curXrxns2set = tmpXrxns2set;
            curXset2 = tmpXset2;
            
            curRxnSet = newRDB;
            curGrowthConditions = newRDB.growthConditions;
            curNonGrowthConditions = newRDB.nonGrowthConditions;
        else
            iterate = 0;
        end
        
    catch gurobiError
        
        if verbose > 1
            gurobiError.message
            fprintf('Error reported during the Expand step\n');
        end
        
        if numIterationsCompleted > 0
            feasible = 0;
        else
            feasible = 1;
        end
        iterate = 0;
    end
end

end

