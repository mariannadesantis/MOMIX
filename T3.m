function [X0,Y0,p,g0,delta,epsilon,epsiloncomp, Q,c,Acon,bcon,Qcon,qcon,ccon] = T3(n) %name of the test instance

%number of objective functions
p=2;
m=2;
r=m+n;

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
c{1} = [1, zeros(1,r-1)];

%second objective function:
const=10;
Q{2} = diag([0,0,const*ones(1,n)]);
c{2} = [0,1,-0.8*const*ones(1,n)];

%linear constraints 
Acon = zeros(1,r);
bcon = 0;

%quadratic constraint
Qcon{1}=diag([ones(1,r)]);
qcon{1}=zeros(r,1);
ccon{1}=4;

%convergence parameters
delta=0.1;
epsilon=0.1;

%parameter for handling numeric rounding errors
epsiloncomp=1e-6;

%box constraints for continous variables
X0=infsup([-2, -2],[2, 2]);
%box constraints for integer variables
Y0=infsup(-2*ones(1,n),2*ones(1,n));



%%%%%%%%%% Functions as they need to be defined for MATLAB %%%%%%%%%%%%%%
%global defined objective functions
global func
func=@(x,i)(1:p==i)*[(x*Q{1})*x' + c{1}*x'; (x*Q{2})*x' + c{2}*x'+n*const*0.16];


%convex constraints
g0=@(x)(x*Qcon{1})*x' + qcon{1}'*x'-ccon{1};

end
