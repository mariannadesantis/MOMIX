function [hp,ListPNS,ListLUB,exitflag,omega_opt2]=integerHyperplanes(lambda,Omega,func,p,m,n,g,ListPNS,ListLUB,epsiloncomp,Q,c,Acon,bcon,Qcon,qcon,ccon)
% function called by Boxcheck2.m in case a further hyperplane to enrich the 
% lower bound considering should be computed, starting from the solution of
% a single-objective mixed integer subproblem.

type=[];
for i=1:m
    type=[type,'C'];
end
for i=m+1:m+n
    type=[type,'I'];
end

% Solve the single-objective mixed integer subproblem in order to compute 
% the further integer hyperplane
[omega_opt2,exitflag]=GurobiCall(lambda,Omega,p,Q,c,Acon,bcon,Qcon,qcon,ccon,type);
add=0;
if exitflag
    % Check if the integer feasible solution found by the solver (omega_opt2) 
    % is an efficient solution:
    hp=f(omega_opt2,func,p);
    [ListPNS, add] = updateListPNS(omega_opt2,m,func,p,g,ListPNS,epsiloncomp);
end

% If omega_opt2 and its function value have been added to ListPNS, 
% the list LUB of local upper bounds is updated
if add
    ListLUB=lub3(ListLUB,f(omega_opt2,func,p),p);
end                

end