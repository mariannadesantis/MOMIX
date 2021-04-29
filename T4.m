function [X0,Y0,p,g0,delta,epsilon,epsiloncomp, Q,c,Acon,bcon,Qcon,qcon,ccon] = T4(m,n) %name of the test instance

%number of objective functions
p=2;
r=2*m+n;

%variable: row vector of length m+n
%Q positive semidefinite (m+n)x(m+n)-matrix
%c row vector of length m+n
%Acon (#constraints)x(m+n) matrix
%bcon column vector of length (#constraints)
%Qcon{i} positive semidefinie (m+n)x(m+n)-matrix
%qcon{i} row vector of length m+n
%ccon{i} scalar

%first objective function:
Q{1} = zeros(r);
c{1} = [ones(1,m), zeros(1,m), ones(1,n)];

%second objective function:
Q{2} = zeros(r);
c{2} = [zeros(1,m), ones(1,m), -ones(1,n)];

%linear constraints 
Acon = zeros(1,r);
bcon = 0;

%quadratic constraint
Qcon{1}=diag([ones(1,2*m), zeros(1,n)]);
qcon{1}=zeros(r,1);
ccon{1}=1;

%convergence parameters
delta=0.1;%0.5;
epsilon=0.1;

%parameter for handling numeric rounding errors
epsiloncomp=1e-3;

%box constraints for continous variables
X0=infsup(-2*ones(1,2*m),2*ones(1,2*m));
%box constraints for integer variables
Y0=infsup(-2*ones(1,n),2*ones(1,n));

%%%%%%%%%% Functions as they need to be defined for MATLAB %%%%%%%%%%%%%%
%global defined objective functions
global func
func=@(x,i)(1:p==i)*[(x*Q{1})*x' + c{1}*x'; (x*Q{2})*x' + c{2}*x'];

%convex constraints
g0=@(x)(x*Qcon{1})*x' + qcon{1}'*x'-ccon{1};

end
