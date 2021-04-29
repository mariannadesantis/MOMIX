function [c,ceq] = P_lubconstr(fun,z,g,lub,p)
% This function is needed to properly call fmincon
omega=z(1:end-1);
t=z(end);
c=[f(omega,fun,p)-lub-t*ones(1,p),g(omega)'];
ceq=[];

end

