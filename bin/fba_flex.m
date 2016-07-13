function [growth,x] = fba_flex(model,exchangeRxnList,growthConditionLB,verbose)
%-------------------------------------------------------------------------- 
% fba_flex - Solves a flux balanc analysis (FBA) problem given a model and 
% a list of lower bounds.
%
% Inputs:
%     model - A metabolic network reconstructions in COBRA format
%     exchangeRxnList - A cell array of exchange reaction IDs which
%       correspond to the rows of the "growthConditionLB" matrix
%     growthConditionLB - A set of lower bounds for the reactions listed in
%       "exchangeRxnList"
%     verbose - 1 indicates that information about the FBA function will be
%       shown
%
% Outputs:
%     growth - The flux value through the biomass rxn
%     x - The flux distribution through the network (one value for each rxn)
%
% Written by Matt Biggs
%--------------------------------------------------------------------------

% Find correspondence between exchange reaction indices
ex_i_m = find(ismember(model.rxns,exchangeRxnList));
ex_i_gcs = zeros(size(ex_i_m));

for k = 1:length(ex_i_m)
    curExRxn = model.rxns{ex_i_m(k)};
    ex_i_gcs(k) = find(ismember(exchangeRxnList,curExRxn));
end

% Set growth conditions
model.lb = zeros(size(model.lb));
model.lb(model.rev > 0) = -1000;
model.lb(ex_i_m) = growthConditionLB(ex_i_gcs);
model.ub = 1000*ones(size(model.lb));

if ~isfield(model,'b')
    model.b = zeros(size(model.mets));
end

% FBA
% Set up primal linear program for Gurobi
sol_primal = struct;
n_mets = size(model.S,1);
n_rxns = size(model.S,2);

A_p = [sparse(model.S); speye(n_rxns); speye(n_rxns)];
RHS_p = [model.b(:); model.ub(:); model.lb(:)];
vtype_p = repmat('C',[n_rxns,1]);
sense_p = [repmat('=',[n_mets,1]); repmat('<',[n_rxns,1]); repmat('>',[n_rxns,1])];
objective_p = model.c;

clear gurobi_model;
gurobi_model.A = A_p;
gurobi_model.obj = objective_p;
gurobi_model.lb = -1e30*ones(size(model.c));
gurobi_model.rhs = RHS_p;
gurobi_model.sense = sense_p;
gurobi_model.vtype = vtype_p;
gurobi_model.modelsense = 'max';
clear params;
params.outputflag = 0;

try
    sol_p = gurobi(gurobi_model, params);
    growth = sol_p.x(model.c > 0);
    x = sol_p.x;
catch gurobiError
    growth = -1;
    x = -1*ones(size(model.c));
    if verbose > 0
        gurobiError
        fprintf('Gurobi error reported in fba_flex\n');
    end
end

end