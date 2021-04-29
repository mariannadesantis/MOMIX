function [X0,Y0,p,g0,delta,epsilon,epsiloncomp,Q,c,A,b,Qcon,qcon,ccon] = T6() %name of the test instance

% not quadratic objective function!
% set all inputs for Gurobi empty -> algorithm2 will recognize it and set
% callMIPsolver=0

%number of objective functions
p=2;
%global defined objective function
global func
func=@(x,i)(1:2==i)*[x(1)+x(3);x(2)+exp(-x(3))];

%box constraints for continous variables
X0=infsup([-2,-2],[2,2]);
%box constraints for integer variables
Y0=infsup([-2],[2]);

%convex constraints
g0=@(x)[x(1)^2+x(2)^2-1];

Q=[];
c=[];
A=[];
b=[];
Qcon=[];
qcon=[];
ccon=[];

%convergence parameters
delta=0.1;
epsilon=0.1;

%parameter for handling numeric rounding errors
epsiloncomp=1e-6;
end
