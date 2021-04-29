function List=sort_list(List)
% This function is needed by MOMIX to sort the boxes 
% within the working list, in order to choose the box with 
% lexicographic smallest ideal point

p=length(List(1).idealpoint);

ind=length(List);
a1=List(end).idealpoint;

for i=1:length(List)-1
    a2=List(i).idealpoint;
    hier=true;
    for j=1:p
        if a2(j)<a1(j)
            hier=false;
            break;
        elseif a1(j)<a2(j)
            hier=true;
            break;
        end
    end
    if hier
        ind=i;
        break;
    end
end
index=[1:ind-1,length(List),ind:(length(List)-1)];
List=List(index);

end

