function [c,ceq] = constr(g,omega,m)
%evaluate the constraint g at point omega
global zul
if zul
    g=@(omega)[];
end
c=g(omega);
ceq=[];
end