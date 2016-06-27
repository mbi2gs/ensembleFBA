function [growth,x] = fba_flex(model,exchangeRxnList,growthConditionLB,verbose)
% Set the lower bounds of model to simulate the input growth condition. The
% lower bounds correspond to the exchange reactions listed in 
% exchangeRxnList. Run FBA and return the flux through the biomass function
% ("growth") and the complete flux distribution ("x").
%
% Written by Matt Biggs, 2016

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