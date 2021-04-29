function [ListPNS, add] = updateListPNS(omega,m, func,p,g, ListPNS, epsiloncomp)
% This function is needed to keep ListPNS (List of Potentially
% Nondominated Solutions) updated along the iterations of MOMIX

r=length(omega);

global feas
%Check, whether omega - given as input - is feasible for MOMINLP:
if all(round(omega(m+1:end))==omega(m+1:end)) &&( feas || all(g(omega)<=0) )
    %If omega is feasible calculate fomega = function value of omega.
    fomega = f(omega,func,p);
    %If fomega is dominated by a function value within ListPNS, 
    % then omega and fomega are not added to ListPNS.
    if ~isempty(ListPNS) && any(all((repmat(fomega,size(ListPNS,1),1)>=ListPNS(:,r+1:end)-epsiloncomp)'))
        add=false;
    else
        %Else, we add omega and fomega to ListPNS
        add=true;
        ListPNS=[ListPNS;omega,fomega];
        %and organize ListPNS in lexicographical order (with respect to the
        %function values):
        [~,index]=sortrows(ListPNS(:,r+1:end));
        ListPNS=ListPNS(index,:);
        %Delete all entries within ListPNS, whose function values are
        %dominated by fomega:
        ListPNS=ListPNS(~min([all((repmat(fomega,size(ListPNS,1),1)<=ListPNS(:,r+1:end))');any((repmat(fomega,size(ListPNS,1),1)<ListPNS(:,r+1:end))')]),:); 
    end
    %If omega is not feasible for MOMINLP, we do not add it to ListPNS.
else
    add=0;
end