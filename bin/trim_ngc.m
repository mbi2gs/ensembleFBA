function [model, rxnDatabase, fromNGCthresh] = phen2net_SEED_trim_ngc(universalRxnSet,growthConditions,nonGrowthConditions,biomassFn,Urxns2set,Uset2,Xrxns2set,Xset2,verbose)
%-------------------------------------------------------------------------- 
% Implementation of the Phen2Net algorithm 
%
% Trim all inconsistent non-growth conditions.
%
% Written by Matt Biggs, mb3ad@virginia.edu, 2016
%-------------------------------------------------------------------------- 
if verbose > 0
    fprintf('\nTrim Non-Growth Condition Step\n');
end

% Initialize outputs
model = struct;
rxnDatabase = struct;
fromNGCthresh = [];

% Append biomass reaction to U
Ubm = [universalRxnSet.S biomassFn(:)];

% Set X (exchange reaction) matrix
X = universalRxnSet.X;
Ex_names = universalRxnSet.Ex_names;

% Identify total number of conditions
n_gc = size(growthConditions,2);
n_ngc = size(nonGrowthConditions,2);
n_cond = n_gc + n_ngc;

% Initialize variables
[n_mets,n_Urxns] = size(Ubm);
[~,n_Xrxns] = size(X);
c_bm = zeros(n_Urxns+n_Xrxns,1); 
c_bm(n_Urxns) = 1;
n_UX = n_Urxns+n_Xrxns;
n_duals = n_mets + n_UX*2; % Includes lambda_mets, labmda_ub, and lambda_lb

% Set proper bounds
ub_U = 1000*ones(n_Urxns,1); % plus one for the biomass rxn
lb_U = [-1000*universalRxnSet.rev; 0];
ub_X = 1000*ones(n_Xrxns,1);

% Build objective
objective = [zeros(1,n_UX*n_cond + n_duals*n_ngc) ones(1,n_ngc)];
vtype = [repmat('C',[1,n_UX*n_cond + n_duals*n_ngc]) repmat('B',[1,n_ngc])];

%## Build matrices
% Build UX diagonal matrix
UX = sparse([Ubm X]);
UX_list = repmat({UX},n_cond,1);
UX_diag = blkdiag(UX_list{:}); % Repeates [Ubm X] along the diagonal n_ngc times

if verbose > 1
    fprintf('Built UX\n');
end

% Build upper and lower bound matrices
lb_mat = sparse(n_UX*n_cond,n_UX*n_cond+n_duals*n_ngc+n_ngc);
lb_rhs = zeros(n_UX*n_cond,1);
for i = 1:n_cond
    if i <= n_gc
        tmp_rhs = [lb_U; growthConditions(:,i)];
    else
        tmp_rhs = [lb_U; nonGrowthConditions(:,i-n_gc)];
    end
    lb_mat(((i-1)*n_UX+1):i*n_UX,:) = [repmat(sparse(n_UX,n_UX),[1,i-1]) speye(n_UX) repmat(sparse(n_UX,n_UX),[1,n_cond-i]) sparse(n_UX,n_duals*n_ngc) sparse(n_UX,n_ngc)];
    lb_rhs(((i-1)*n_UX+1):i*n_UX,1) = tmp_rhs;
end

ub_mat = sparse(n_UX*n_cond,n_UX*n_cond+n_duals*n_ngc+n_ngc);
ub_rhs = repmat([ub_U; ub_X],n_cond,1);
for i = 1:n_cond
    ub_mat(((i-1)*n_UX+1):i*n_UX,:) = [repmat(sparse(n_UX,n_UX),[1,i-1]) speye(n_UX) repmat(sparse(n_UX,n_UX),[1,n_cond-i]) sparse(n_UX,n_duals*n_ngc) sparse(n_UX,n_ngc)];
end

if verbose > 1
    fprintf('Built LB/UB mats\n');
end

% Build dual problem constraint matrices
dualc_mat = sparse(n_UX*n_ngc,n_UX*n_cond+n_duals*n_ngc+n_ngc);
for i = 1:n_ngc
    dualc_mat(((i-1)*n_UX+1):i*n_UX,:) = [sparse(n_UX,n_UX*n_cond) repmat(sparse(n_UX,n_duals),[1,i-1]) UX' speye(n_UX) -speye(n_UX) repmat(sparse(n_UX,n_duals),[1,n_ngc-i]) sparse(n_UX,n_ngc)];
end

dual_sign_mat = sparse(n_UX*n_ngc*2,n_UX*n_cond+n_duals*n_ngc+n_ngc);
for i = 1:n_ngc
    dual_sign_mat(((i-1)*n_UX*2+1):i*n_UX*2,:) = [sparse(n_UX*2,n_UX*n_cond) repmat(sparse(n_UX*2,n_duals),[1,i-1]) sparse(n_UX*2,n_mets) speye(n_UX*2) repmat(sparse(n_UX*2,n_duals),[1,n_ngc-i]) sparse(n_UX*2,n_ngc)];
end

dual_obj_mat = sparse(n_ngc,n_UX*n_cond+n_duals*n_ngc+n_ngc);
for i = 1:n_ngc
    tmp_ub = sparse([ub_U(:)' ub_X(:)']);
    tmp_lb = -sparse(double([lb_U(:)' nonGrowthConditions(:,i)']));
    dual_obj_mat(i,:) = [sparse(1,n_UX*n_gc) repmat(sparse(1,n_UX),[1,i-1]) -c_bm(:)' repmat(sparse(1,n_UX),[1,n_ngc-i]) repmat(sparse(1,n_duals),[1,i-1]) sparse(1,n_mets) tmp_ub tmp_lb repmat(sparse(1,n_duals),[1,n_ngc-i]) sparse(1,n_ngc)];
end

if verbose > 1
    fprintf('Built dual mats\n');
end

% Build biomass threshold matrices
c_gc_list = repmat({c_bm'},n_gc,1);
c_gc_diag = sparse(blkdiag(c_gc_list{:}));
c_ngc_list = repmat({c_bm'},n_ngc,1);
c_ngc_diag = sparse(blkdiag(c_ngc_list{:}));

% Put it all together
A =  [UX_diag sparse(n_mets*n_cond,n_duals*n_ngc+n_ngc) ; ...                              %(1) Steady state assumption     
      lb_mat; ...                                                                               %(2) Flux bounds
      ub_mat; ...                                                                               % 
      c_gc_diag sparse(n_gc,n_UX*n_ngc+n_duals*n_ngc+n_ngc); ...                           %(3) Constrain gc biomass > 0.05 
      sparse(n_ngc,n_UX*n_gc) c_ngc_diag sparse(n_ngc,n_duals*n_ngc) 1000*speye(n_ngc);... %(4) Upper bound on ngc biomass flux = 0 unless that ngc is removed
      dualc_mat; ...                                                                            %(5) Duality variable constraints      
      dual_sign_mat; ...                                                                        %    Dual variables are strictly positive (except lambda_mets)
      dual_obj_mat; ...                                                                         %(6) require that primal and dual objectives equal each other
      ]; 

clear UX_diag lb_mat ub_mat lb_mat_z ub_mat_z lb_z c_gc_diag c_ngc_diag c_gc_list c_ngc_list dualc_mat UX UX_list dual_sign_mat y_relax_rs_sign_mat y_relax_lambda_lb_sign_mat

if verbose > 1
    fprintf('Built A\n');
end

RHS = [ zeros(n_mets*n_cond,1); ...                     %(1)  
        lb_rhs; ...                                     %(2)
        ub_rhs; ...                                     %
        0.05*ones(n_gc,1); ...                          %(3)
        1000*ones(n_ngc,1); ...                         %(4)
        repmat(c_bm,n_ngc,1); ...                       %(5)        
        zeros(n_UX*n_ngc*2,1); ...                      %
        zeros(n_ngc,1); ...                             %(6)
        ];                           

sense = [repmat('=',n_mets*n_cond,1); ...               %(1)    
         repmat('>',n_UX*n_cond,1); ...                 %(2)
         repmat('<',n_UX*n_cond,1); ...                 %
         repmat('>',n_gc,1); ...                        %(3)
         repmat('<',n_ngc,1); ...                       %(4)
         repmat('=',n_UX*n_ngc,1);...                   %(5)         
         repmat('>',n_UX*n_ngc*2,1);...                 %
         repmat('=',n_ngc,1); ...                       %(6)
         ];                                             

if verbose > 1     
    fprintf('Built entire problem\n');     
end

clear dual_sign_mat dual_sign_mat_z lb_rhs

if verbose > 1
    fprintf(['\nObjective size:   ' num2str(size(objective)) '\n']);
    fprintf(['Vtype size:       ' num2str(size(vtype)) '\n']);
    fprintf(['A size:           ' num2str(size(A)) '\n']);
    fprintf(['RHS size:         ' num2str(size(RHS)) '\n']);
    fprintf(['Sense size:       ' num2str(size(sense)) '\n\n']);
end

try
    clear gurobi_model;
    gurobi_model.A = A;
    gurobi_model.obj = objective;
    gurobi_model.lb = [-1e30*ones(1,n_UX*n_cond + n_duals*n_ngc) zeros(1,n_ngc)];
    gurobi_model.rhs = RHS;
    gurobi_model.sense = sense;
    gurobi_model.vtype = vtype;
    gurobi_model.modelsense = 'max';

    clear params;
    params.outputflag = 0;

    result = gurobi(gurobi_model, params);

    if verbose > 1
        disp(result)
    end

    n_front = n_UX*n_cond + n_duals*n_ngc;
    r_ngc = result.x(n_front+1 : n_front+n_ngc);
    rngclog = r_ngc > 1e-10;

    % Make COBRA model of current solution
    newModel = struct;
    newModel.mets = universalRxnSet.mets;
    newModel.metFormulas = universalRxnSet.metFormulas;
    newModel.metNames = universalRxnSet.metNames;
    newModel.rxns = [universalRxnSet.rxns; Ex_names; {'Biomass'}];
    newModel.rxnNames = [universalRxnSet.rxnNames; Ex_names; {'Biomass'}];
    newModel.rev = double([universalRxnSet.rev; zeros(length(Ex_names),1); 0]);
    newModel.S = double([universalRxnSet.S X biomassFn(:)]);
    newModel.c = double(zeros(n_UX,1)); newModel.c(end,1) = 1;
    newModel.b = double(zeros(length(universalRxnSet.mets),1));
    newModel.ub = double([ub_U; ub_X; 1000]);
    newModel.lb = double([lb_U; growthConditions(:,1); 0]);
    model = newModel;

    % Format as a reaction database for the expansion algorithm        
    newRDB = struct;
    newRDB.mets = universalRxnSet.mets;
    newRDB.metFormulas = universalRxnSet.metFormulas;
    newRDB.metNames = universalRxnSet.metNames;
    newRDB.rxns = universalRxnSet.rxns;
    newRDB.rxnNames = universalRxnSet.rxnNames;
    newRDB.rev = universalRxnSet.rev;
    newRDB.S = universalRxnSet.S;
    newRDB.X = X;
    newRDB.Ex_names = Ex_names;
    newRDB.growthConditions = growthConditions;
    newRDB.nonGrowthConditions = nonGrowthConditions(:,rngclog);
    rxnDatabase = newRDB;

    fromNGCthresh = rngclog;
    
catch gurobiError
    gurobiError.message
    fprintf('Error reported\n');
end

end

