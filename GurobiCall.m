function [omega_opt2,exitflag]=GurobiCall(lambda,Omega,p,Q,c,Acon,bcon,Qcon, qcon, ccon, type)

global GurobiTime
 
l=inf(Omega); 
u=sup(Omega); 

clear model;
Qobj = lambda(1)*Q{1};
cobj = lambda(1)*c{1};
for j=2:p
    Qobj=Qobj+lambda(j)*Q{j};
    cobj=cobj+lambda(j)*c{j};
end

%Add objective function
model.Q = sparse(Qobj);
model.obj = cobj;


%linear constraints
model.A = sparse(Acon);
model.rhs = bcon;
model.sense ='<';


% quadratic constraints
for j=1:length(Qcon)
    model.quadcon(j).Qc=sparse(Qcon{j});
    model.quadcon(j).q=qcon{j};
    model.quadcon(j).rhs=ccon{j};
    model.quadcon(j).sense='<';
    model.quadcon(j).name=['quadratic_constraint_',num2str(j)];
end

model.ub = u;
model.lb = l;

% Integrality constraints
model.vtype =type;
model.modelsense = 'min';
% model.varnames = names;

gurobi_write(model, 'miqp1.lp');

clear params;
params.outputflag = 0;
gurobistart=toc;
result = gurobi(model, params);
GurobiTime=GurobiTime+toc-gurobistart;
%result.status
omega_opt2=1:size(Omega,2);

if strcmp(result.status,'OPTIMAL')
    exitflag = 1;
    for v=1:size(Omega,2)
        omega_opt2(v) = result.x(v);
    end
elseif strcmp(result.status,'INFEASIBLE')
    exitflag = -2;
    return
else
    exitflag = -1;
    for v=1:size(Omega,2)
        omega_opt2(v) = result.x(v);
    end    
end

end
    
    
