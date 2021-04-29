function z = f(x,func,p)
% evaluate the objective function at point x:
for i=1:p
    z(i)=func(x,i);
end
end
