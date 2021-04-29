function ListLUBnew = lub3(ListLUB,z,p)
%This function is the procedure to update the local upper bound set 
%from the list of potentially nondominated solutions

A=[];ListLUBnew=[];
for i=1:size(ListLUB,1)
    u=ListLUB(i,:);
    if all(z<u)
        A=[A;u];
    end
end
B=cell(p,1);P=B;Ph=P;

for j=1:p
    Bj=B{j};
    for i=1:size(ListLUB,1)
        u=ListLUB(i,:);
        if z(j)==u(j) && all([z(1:j-1),z(j+1:p)]<[u(1:j-1),u(j+1:p)])
            Bj=[Bj;u];
        end
    end
    B{j}=Bj;
end

for i=1:size(A,1)
    u=A(i,:);
    for j=1:p
        Pj=P{j};
        Pj=[Pj;u(1:j-1),z(j),u(j+1:p)];
        P{j}=Pj;
    end
end
V=cell(p,1);
for j=1:p
    V{j}=[P{j};B{j}];
end
for j=1:p
    Vj=V{j};Pj=P{j};Bj=B{j};Phj=Ph{j};
    for i=1:size(Pj,1)
        y=Pj(i,:);
        redundant=false;
        for k=1:size(Vj,1)
            u=Vj(k,:);
            if all(y<=u) && any(y<u)
                redundant=true;break;
            end
        end
        if ~redundant
            Phj=[Phj;y];
        end
    end
    V{j}=Vj;P{j}=Pj;B{j}=Bj;Ph{j}=Phj;
end

for i=1:size(ListLUB,1)
    u=ListLUB(i,:);
    if isempty(A)||isempty(finde( A, u)) 
        ListLUBnew=[ListLUBnew;u];
    end
end
for j=1:p
    Phj=Ph{j};
    for k=1:size(Phj,1)
        u=Phj(k,:);
        if isempty(finde(ListLUBnew,u))
            ListLUBnew=[ListLUBnew;u];
        end
    end
end