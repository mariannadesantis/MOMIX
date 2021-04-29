function [X0,Y0,p,g0,delta,epsilon,epsiloncomp,Q,c,A,b,Qcon,qcon,ccon] = T5() %name of the test instance

%number of objective functions
p=3;

%variable: row vector of length m+n
%Q positive semidefinite (m+n)x(m+n)-matrix
%c row vector of length m+n
%Acon (#constraints)x(m+n) matrix
%bcon column vector of length (#constraints)
%Qcon{i} positive semidefinie (m+n)x(m+n)-matrix
%qcon{i} row vector of length m+n
%ccon{i} scalar

%first objective function:
Q{1} = zeros(4);
c{1} = [1,0,0,1];

%second objective function:
Q{2} = zeros(4);
c{2} = [0,1,0,-1];

%third objective function:
Q{3} = diag([0,0,0,1]);
c{3} = [0,0,1,0];

%linear constraints 
A = zeros(1,4);
b = 0;

%quadratic constraint
Qcon{1}=diag([1,1,1,0]);
qcon{1}=zeros(4,1);
ccon{1}=1;

%global defined objective function
global func
func=@(x,i)(1:p==i)*[(x*Q{1})*x' + c{1}*x'; (x*Q{2})*x' + c{2}*x'; (x*Q{3})*x' + c{3}*x'];

%box constraints for continous variables
X0=infsup([-2,-2,-2],[2,2,2]);
%box constraints for integer variables
Y0=infsup([-2],[2]);

%convex constraints
g0=@(x) (x*Qcon{1})*x' + qcon{1}'*x'-ccon{1};


%convergence parameters
delta=0.5;
epsilon=0.1;

%parameter for handling numeric rounding errors
epsiloncomp=1e-6;
end
