function [model, rxnDatabase, fromDecthresh, fromGCthresh] = trim_active(rxnSet,growthConditions,nonGrowthConditions,biomassFn,decisionRxns,verbose)
%-------------------------------------------------------------------------- 
% trim_active - Removes reactions from rxnSet in order to prevent growth in
% the nonGrowthConditions.
%
% Inputs:
%     rxnSet - Matlab structure containing an S matrix and a similar
%       matrix for exchange rxns (X matrix), a reversability indicator for
%       all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange 
%       rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames), 
%       and metabolite formulas (metFormulas).
%     growthConditions - set of lower bounds corresponding to growth media conditions
%     nonGrowthConditions - set of lower bounds corresponding to non-growth media conditions
%     biomassFn - same format as a reaction in universalRxnSet.S
%     decisionRxns - List of rxn indices (from rxnSet.U) that will be candidates for trimming
%     verbose - Indicates verbose output (1 is some detail, 2 is all details)
%
% Outputs:
%     model - A metabolic network reconstruction in COBRA format as a matlab struct
%     rxnDatabase - The model in the same format as universalRxnSet
%     fromDecthresh - logical vector indicating the rxns kept (1) and removed (0)
%     fromGCthresh - logical vector indicating the growth conditions kept (1) and removed (0)
%
% Written by Matt Biggs, mb3ad@virginia.edu, 2016
%-------------------------------------------------------------------------- 
if verbose > 0
    fprintf('\nTrim Step\n');
end

% Initialize outputs
model = struct;
rxnDatabase = struct;
fromDecthresh = [];
fromGCthresh = [];

% Append biomass reaction to U
Ubm = [rxnSet.S biomassFn(:)];

% Set X (exchange reaction) matrix
X = rxnSet.X;
Ex_names = rxnSet.Ex_names;

% Identify total number of conditions
n_gc = size(growthConditions,2);
n_ngc = size(nonGrowthConditions,2);
n_cond = n_gc + n_ngc;
n_decs = length(decisionRxns);
gcRemovalCost = 100;

% Initialize variables
[n_mets,n_Urxns] = size(Ubm);
[~,n_Xrxns] = size(X);
c_bm = zeros(n_Urxns+n_Xrxns,1); 
c_bm(n_Urxns) = 1;
n_UX = n_Urxns+n_Xrxns;
n_duals = n_mets + n_UX*2 + n_decs*2; % Includes lambda_mets, labmda_ub, yrelax_ub, lambda_lb and yrelax_lb
dec_log = zeros(n_UX,1);
dec_log(decisionRxns) = 1; 
dec_log = dec_log > 0;

% Set proper bounds
ub_U = 1000*ones(n_Urxns,1); % plus one for the biomass rxn
lb_U = [-1000*rxnSet.rev; 0];
ub_X = 1000*ones(n_Xrxns,1);

% Build objective
objective = [zeros(1,n_UX*n_cond + n_duals*n_ngc) ones(1,n_decs) gcRemovalCost*ones(1,n_gc)];
vtype = [repmat('C',[1,n_UX*n_cond + n_duals*n_ngc]) repmat('B',[1,n_decs+n_gc])];

%## Build matrices
% Build UX diagonal matrix
UX = sparse([Ubm X]);
UX_list = repmat({UX},n_cond,1);
UX_diag = blkdiag(UX_list{:}); % Repeates [Ubm X] along the diagonal n_ngc times

if verbose > 1
    fprintf('Built UX\n');
end

% Build upper and lower bound matrices
lb_dec_mat = sparse(n_UX*n_cond,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
lb_dec_rhs = zeros(n_UX*n_cond,1);
for i = 1:n_cond    
    if i <= n_gc
        tmp_rhs = [lb_U; growthConditions(:,i)];
    else
        tmp_rhs = [lb_U; nonGrowthConditions(:,i-n_gc)];
    end
    
    tmp_ind_mat = sparse(n_UX,n_decs);
    indices = sub2ind(size(tmp_ind_mat),decisionRxns(:)',1:length(decisionRxns));
    tmp_ind_mat(indices) = -tmp_rhs(dec_log);
    tmp_rhs(dec_log) = 0;
     
    lb_dec_mat(((i-1)*n_UX+1):i*n_UX,:) = [repmat(sparse(n_UX,n_UX),[1,i-1]) speye(n_UX) repmat(sparse(n_UX,n_UX),[1,n_cond-i]) sparse(n_UX,n_duals*n_ngc) tmp_ind_mat sparse(n_UX,n_gc)];
    lb_dec_rhs(((i-1)*n_UX+1):i*n_UX,1) = tmp_rhs;
end

ub_dec_mat = sparse(n_UX*n_cond,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
ub_dec_rhs = zeros(n_UX*n_cond,1);
for i = 1:n_cond
    tmp_rhs = [ub_U(:);ub_X(:)];
    
    tmp_ind_mat = sparse(n_UX,n_decs);
    indices = sub2ind(size(tmp_ind_mat),decisionRxns(:)',1:length(decisionRxns));
    tmp_ind_mat(indices) = -tmp_rhs(dec_log);
    tmp_rhs(dec_log) = 0;
    
    ub_dec_mat(((i-1)*n_UX+1):i*n_UX,:) = [repmat(sparse(n_UX,n_UX),[1,i-1]) speye(n_UX) repmat(sparse(n_UX,n_UX),[1,n_cond-i]) sparse(n_UX,n_duals*n_ngc) tmp_ind_mat sparse(n_UX,n_gc)];
    ub_dec_rhs(((i-1)*n_UX+1):i*n_UX,1) = tmp_rhs;
end

if verbose > 1
    fprintf('Built LB/UB mats\n');
end

% Build dual problem constraint matrices
dualc_mat = sparse(n_UX*n_ngc,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
for i = 1:n_ngc
    dualc_mat(((i-1)*n_UX+1):i*n_UX,:) = [sparse(n_UX,n_UX*n_cond) repmat(sparse(n_UX,n_duals),[1,i-1]) UX' speye(n_UX) sparse(n_UX,n_decs) -speye(n_UX) sparse(n_UX,n_decs) repmat(sparse(n_UX,n_duals),[1,n_ngc-i]) sparse(n_UX,n_decs) sparse(n_UX,n_gc)];
end

rowCnt = (n_UX*2+n_decs*2);
dual_sign_mat = sparse(rowCnt*n_ngc,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
for i = 1:n_ngc
    dual_sign_mat(((i-1)*rowCnt+1):i*rowCnt,:) = [sparse(rowCnt,n_UX*n_cond) repmat(sparse(rowCnt,n_duals),[1,i-1]) sparse(rowCnt,n_mets) speye(rowCnt) repmat(sparse(rowCnt,n_duals),[1,n_ngc-i]) sparse(rowCnt,n_decs) sparse(rowCnt,n_gc)];
end

Uub = 1000;
y_relax_rs_ub_mat = sparse(n_decs*n_ngc,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
for i = 1:n_ngc    
    y_relax_rs_ub_mat(((i-1)*n_decs+1):i*n_decs,:) = [sparse(n_decs,n_UX*n_cond) repmat(sparse(n_decs,n_duals),[1,i-1]) sparse(n_decs,n_mets+n_UX) speye(n_decs) sparse(n_decs,n_UX+n_decs) repmat(sparse(n_decs,n_duals),[1,n_ngc-i]) -Uub*speye(n_decs) sparse(n_decs,n_gc)];
end

y_relax_rs_lb_mat = sparse(n_decs*n_ngc,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
for i = 1:n_ngc
    y_relax_rs_lb_mat(((i-1)*n_decs+1):i*n_decs,:) = [sparse(n_decs,n_UX*n_cond) repmat(sparse(n_decs,n_duals),[1,i-1]) sparse(n_decs,n_mets+n_UX*2+n_decs) speye(n_decs) repmat(sparse(n_decs,n_duals),[1,n_ngc-i]) -Uub*speye(n_decs) sparse(n_decs,n_gc)];
end

y_relax_lambda_ub_sign_mat = sparse(n_decs*n_ngc,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
for i = 1:n_ngc
    tmp_ub = -spdiags([ub_U(:); ub_X(:)],0,n_UX,n_UX);
    y_relax_lambda_ub_sign_mat(((i-1)*n_decs+1):i*n_decs,:) = [sparse(n_decs,n_UX*n_cond) repmat(sparse(n_decs,n_duals),[1,i-1]) sparse(n_decs,n_mets) tmp_ub(dec_log,:) speye(n_decs) sparse(n_decs,n_UX+n_decs) repmat(sparse(n_decs,n_duals),[1,n_ngc-i]) sparse(n_decs,n_decs) sparse(n_decs,n_gc)];
end

y_relax_lambda_lb_sign_mat = sparse(n_decs*n_ngc,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
for i = 1:n_ngc
    tmp_lb = spdiags([lb_U(:); nonGrowthConditions(:,i)],0,n_UX,n_UX);
    y_relax_lambda_lb_sign_mat(((i-1)*n_decs+1):i*n_decs,:) = [sparse(n_decs,n_UX*n_cond) repmat(sparse(n_decs,n_duals),[1,i-1]) sparse(n_decs,n_mets) sparse(n_decs,n_UX+n_decs) tmp_lb(dec_log,:) speye(n_decs) repmat(sparse(n_decs,n_duals),[1,n_ngc-i]) sparse(n_decs,n_decs) sparse(n_decs,n_gc)];
end

y_relax_lambda_ub_rs_sign_mat = sparse(n_decs*n_ngc,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
for i = 1:n_ngc
    tmp_ub = -spdiags([ub_U(:); ub_X(:)],0,n_UX,n_UX);
    y_relax_lambda_ub_rs_sign_mat(((i-1)*n_decs+1):i*n_decs,:) = [sparse(n_decs,n_UX*n_cond) repmat(sparse(n_decs,n_duals),[1,i-1]) sparse(n_decs,n_mets) tmp_ub(dec_log,:) speye(n_decs) sparse(n_decs,n_UX+n_decs) repmat(sparse(n_decs,n_duals),[1,n_ngc-i]) -Uub*speye(n_decs) sparse(n_decs,n_gc)];
end

y_relax_lambda_lb_rs_sign_mat = sparse(n_decs*n_ngc,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
for i = 1:n_ngc
     tmp_lb = spdiags([lb_U(:); nonGrowthConditions(:,i)],0,n_UX,n_UX);
    y_relax_lambda_lb_rs_sign_mat(((i-1)*n_decs+1):i*n_decs,:) = [sparse(n_decs,n_UX*n_cond) repmat(sparse(n_decs,n_duals),[1,i-1]) sparse(n_decs,n_mets) sparse(n_decs,n_UX+n_decs) tmp_lb(dec_log,:) speye(n_decs) repmat(sparse(n_decs,n_duals),[1,n_ngc-i]) -Uub*speye(n_decs) sparse(n_decs,n_gc)];
end

dual_obj_mat = sparse(n_ngc,n_UX*n_cond+n_duals*n_ngc+n_decs+n_gc);
for i = 1:n_ngc
    tmp_ub = sparse([ub_U(:); ub_X(:)]);
    tmp_ub(dec_log) = 0;
    tmp_lb = -sparse(double([lb_U(:); nonGrowthConditions(:,i)]));
    tmp_lb(dec_log) = 0;
    dual_obj_mat(i,:) = [sparse(1,n_UX*n_gc) repmat(sparse(1,n_UX),[1,i-1]) -c_bm(:)' repmat(sparse(1,n_UX),[1,n_ngc-i]) repmat(sparse(1,n_duals),[1,i-1]) sparse(1,n_mets) tmp_ub(:)' sparse(ones(1,n_decs)) tmp_lb(:)' sparse(ones(1,n_decs)) repmat(sparse(1,n_duals),[1,n_ngc-i]) sparse(1,n_decs+n_gc)];
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
A =  [UX_diag sparse(n_mets*n_cond,n_duals*n_ngc+n_decs+n_gc) ; ...                 %(1) Steady state assumption     
      lb_dec_mat; ...                                                               %(2) Flux bounds
      ub_dec_mat; ...                                                               % 
      c_gc_diag sparse(n_gc,n_UX*n_ngc+n_duals*n_ngc+n_decs) -0.05*speye(n_gc); ... %(3) Constrain gc biomass > 0.05 
      sparse(n_ngc,n_UX*n_gc) c_ngc_diag sparse(n_ngc,n_duals*n_ngc+n_decs+n_gc);...%(4) Upper bound on ngc biomass flux = 0
      dualc_mat; ...                                                                %(5) Duality variable constraints      
      dual_sign_mat; ...                                                            %    Dual variables are strictly positive (except lambda_mets)
	  y_relax_rs_ub_mat; ...                                                        %(6) Relax quadtratic binary constraints to linear constraints            
      y_relax_rs_lb_mat; ...                                                        %
      y_relax_lambda_ub_sign_mat; ...                                               %
	  y_relax_lambda_lb_sign_mat; ...                                               %
      y_relax_lambda_ub_rs_sign_mat; ...                                            %
	  y_relax_lambda_lb_rs_sign_mat; ...                                            %
      dual_obj_mat; ...                                                             %(7) require that primal and dual objectives equal each other
      ]; 

clear UX_diag lb_mat ub_mat lb_mat_z ub_mat_z lb_z c_gc_diag c_ngc_diag c_gc_list c_ngc_list dualc_mat UX UX_list dual_sign_mat y_relax_rs_sign_mat y_relax_lambda_lb_sign_mat

if verbose > 1
    fprintf('Built A\n');
end

RHS = [ zeros(n_mets*n_cond,1); ...                     %(1)  
        lb_dec_rhs; ...                                 %(2)
        ub_dec_rhs; ...                                 %
        zeros(n_gc,1); ...                              %(3)
        zeros(n_ngc,1); ...                             %(4)
        repmat(c_bm,n_ngc,1); ...                       %(5)
        zeros((n_UX*2+n_decs*2)*n_ngc,1); ...           %
		zeros(n_decs*n_ngc,1); ...                      %(6)
        zeros(n_decs*n_ngc,1); ...                      %
        zeros(n_decs*n_ngc,1); ...                      %
		zeros(n_decs*n_ngc,1); ...                      %
        -Uub*ones(n_decs*n_ngc,1); ...                  %
        -Uub*ones(n_decs*n_ngc,1); ...                  %
        zeros(n_ngc,1); ...                             %(7)
        ];                           

sense = [repmat('=',n_mets*n_cond,1); ...               %(1)    
         repmat('>',n_UX*n_cond,1); ...                 %(2)
         repmat('<',n_UX*n_cond,1); ...                 %
         repmat('>',n_gc,1); ...                        %(3)
         repmat('=',n_ngc,1); ...                       %(4)
         repmat('=',n_UX*n_ngc,1);...                   %(5)
         repmat('>',(n_UX*2+n_decs*2)*n_ngc,1);...      %
		 repmat('<',n_decs*n_ngc,1);...                 %(6)
         repmat('<',n_decs*n_ngc,1);...                 %
         repmat('<',n_decs*n_ngc,1);...                 %         
		 repmat('<',n_decs*n_ngc,1);...                 %
         repmat('>',n_decs*n_ngc,1);...                 %
		 repmat('>',n_decs*n_ngc,1);...                 %
         repmat('=',n_ngc,1); ...                       %(7)
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
    gurobi_model.lb = [-1e30*ones(1,n_UX*n_cond + n_duals*n_ngc) zeros(1,n_decs+n_gc)];
    gurobi_model.rhs = RHS;
    gurobi_model.sense = sense;
    gurobi_model.vtype = vtype;
    gurobi_model.modelsense = 'max';

    clear params;
    params.outputflag = 0;
    params.SolutionLimit = 5;
    params.TimeLimit = 200;
    
    result = gurobi(gurobi_model, params);

    if verbose > 1
        disp(result)
    end
    
    n_front = n_UX*n_cond + n_duals*n_ngc;
    r_decs = result.x(n_front+1 : n_front+n_decs);
    rdecslog = r_decs > 1e-10;
    rulog = ones(n_Urxns-1,1);
    rulog(decisionRxns) = rdecslog;
    rulog = rulog > 0;
    r_gcs = result.x(n_front+n_decs+1 : n_front+n_decs+n_gc);
    rgclog = r_gcs > 1e-10;
    
    % Make COBRA model of current solution
    newModel = struct;
    newModel.mets = rxnSet.mets;
    newModel.metFormulas = rxnSet.metFormulas;
    newModel.metNames = rxnSet.metNames;
    newModel.rxns = [rxnSet.rxns(rulog); Ex_names; {'Biomass'}];
    newModel.rxnNames = [rxnSet.rxnNames(rulog); Ex_names; {'Biomass'}];
    newModel.rev = double([rxnSet.rev(rulog); zeros(length(Ex_names),1); 0]);
    newModel.S = double([rxnSet.S(:,rulog) X biomassFn(:)]);
    newModel.c = double(zeros(sum(rulog) + 1 + length(Ex_names),1)); newModel.c(end,1) = 1;
    newModel.b = double(zeros(length(rxnSet.mets),1));
    newModel.ub = double([ub_U([rulog; 1>10]); ub_X; 1000]);
    newModel.lb = double([lb_U([rulog; 1>10]); growthConditions(:,1); 0]);
    model = newModel;

    % Format as a reaction database for the expansion algorithm        
    newRDB = struct;
    newRDB.mets = rxnSet.mets;
    newRDB.metFormulas = rxnSet.metFormulas;
    newRDB.metNames = rxnSet.metNames;
    newRDB.rxns = rxnSet.rxns(rulog);
    newRDB.rxnNames = rxnSet.rxnNames(rulog);
    newRDB.rev = rxnSet.rev(rulog);
    newRDB.S = rxnSet.S(:,rulog);
    newRDB.X = X;
    newRDB.Ex_names = Ex_names;
    newRDB.growthConditions = growthConditions;
    newRDB.nonGrowthConditions = nonGrowthConditions;
    rxnDatabase = newRDB;

    fromDecthresh = rulog;
    fromGCthresh = rgclog;
    
catch gurobiError
    if verbose > 0
        gurobiError.message
        fprintf('Error reported\n');
    end
end

end

