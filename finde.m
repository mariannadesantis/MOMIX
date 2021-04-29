function ind = finde(A,b)

if isempty(A) || isempty(b)
    ind=[];
    flag=false;
else
    [a,ind]=max(all((repmat(b,size(A,1),1)==A)'));
    flag=true;
end

if flag && ~a
    ind=[];
end

end