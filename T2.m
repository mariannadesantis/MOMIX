function [X0,Y0,p,g0,delta,epsilon,epsiloncomp,Q,c,Acon,bcon,Qcon,qcon,ccon] = T2(r) %name of the test instance

% example file for biulding a quadratic convex test instance which can be
% solved by MOMIX
% Input: r is the  total number of variables

% number of objective functions
p=2;

% variable: row vector of length m+n
% Q positive semidefinite (m+n)x(m+n)-matrix
% c row vector of length m+n
% Acon (#constraints)x(m+n) matrix
% bcon column vector of length (#constraints)
% Qcon{i} positive semidefinie (m+n)x(m+n)-matrix
% qcon{i} row vector of length m+n
% ccon{i} scalar

% first objective function:
Q1 = ones(r)+diag([2,zeros(1,r-2),3]);
% Q{1}=1/r*(Q1'*Q1);
Q{1}=Q1'*Q1;
c{1} = [1,2*ones(1,r-2),1];

% second objective function:
Q2 = ones(r)+diag([1,3*ones(1,r-2),1]);% Also interesting A2=[2 1 1;1 4 1; 1 1 1];
% Q{2}=1/r*(Q2'*Q2);
Q{2}=Q2'*Q2;
c{2} = [-1,-2*ones(1,r-2),5];

% add more objectives by using Q{j}, c{j} for j=3,4,...

% linear constraints 
Acon = zeros(1,r);
bcon = 0;

% quadratic constraint; add more in Qcon{2},....
Qcon{1}=zeros(r);
qcon{1}=zeros(r,1);
ccon{1}=0;

% add more quadratic constraints by using Qcon{j}, qcon{j}, ccon{j}, for
% j=2,3,...

% global defined objective functions - in case, add more objective
% functions
global func
func=@(x,i)(1:p==i)*[(x*Q{1})*x' + c{1}*x'; (x*Q{2})*x' + c{2}*x'];

% box constraints for continous variables
X0=infsup([-5,-5],[5,5]);
% box constraints for integer variables
Y0=infsup(-5*ones(1,r-2),5*ones(1,r-2));
% if purely integer/continuous, use Y0/X0 = infsup([], []);

% convex constraints - in case, add more quadratic contraints; if one kind
% of constraints not existent, remove these
g0=@(x)[];

% convergence parameters
delta=0.1;      % maximal boxwidth of the boxes in the solution list
epsilon=0.1;    %not improtant at the moment

%parameter for handling numeric rounding errors
epsiloncomp=1e-3;
end
