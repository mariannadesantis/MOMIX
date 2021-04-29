function h=plot_image(p,FV,ListPNS)
% Plot the solution in the image/criterion space
if p==2
    if nargout
        h=figure;
    end
    if ~isempty(FV)
        plot(FV(:,end-1),FV(:,end),'LineStyle','none','Marker','.','MarkerEdgeColor',0.7*[1 1 1])
        hold on;
    end    
    plot(ListPNS(:,end-1),ListPNS(:,end),'LineStyle','none','Marker','.','MarkerEdgeColor',[0 0 0])
    hold on;
    xlabel('f1(x)');
    ylabel('f2(x)');
    set(gca,'Color',[0.9 0.9 0.9], 'Fontsize', 15)
elseif p==3
    if nargout
        h=figure;
    end
    if ~isempty(FV)
        scatter3(FV(:,end-2),FV(:,end-1), FV(:,end),'.','MarkerEdgeColor', 0.5*[1 1 1]);
        hold on;
    end
    scatter3(ListPNS(:,end-2),ListPNS(:,end-1),ListPNS(:,end),'.','MarkerEdgeColor',[0 0 0])
    hold on;
    xlabel('f1(x)');
    ylabel('f2(x)');
    zlabel('f3(x)');
    set(gca,'Fontsize', 15)        
else
    h=[];
end

drawnow;
hold off;


end