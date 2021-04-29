function h=plot_boxes(Omega0,ListS,ListNS,ListW)
% Plot the solution in the pre-image/decision space

if length(Omega0)~=2
    if nargout
        h=[];
    end
    return;
else
    for i = 1:length(ListNS)
        flag=ListNS(i).flag;
        switch flag
            case 0 %infeasible box
                ListNS(i).Farbe=0.6*[1 1 1];
            case 1 %discarded by discarding test
                ListNS(i).Farbe=0.7*[1 1 1];
            case 2 %discarded in the second while-loop
                ListNS(i).Farbe=0.8*[1 1 1];
            case 3
                ListNS(i).Farbe='r';
        end
    end
    
    
    ListeS_Farbe=0.3*[1 1 1];
    ListeW_Farbe='y';
        
    a=inf(Omega0(1));
    b=sup(Omega0(1));
    c=inf(Omega0(2));
    d=sup(Omega0(2));    
    
    if nargout
         h=plot([a,a,b,b,a],[c,d,d,c,c], 'k'); 
         hold on;
    end    
    plot([a,a,b,b,a],[c,d,d,c,c], 'k')  
    hold on;
    fill([a,a,b,b,a], [c,d,d,c,c], 'w');
    xlim([inf(Omega0(1)) sup(Omega0(1))]) 
    ylim([inf(Omega0(2)) sup(Omega0(2))]) 
    set(gca,'Color',[0.9 0.9 0.9]) 
    xlabel('x1')
    ylabel('x2')
    set(gca,'FontSize',14)
    
    for i=1:length(ListW)
        Box=ListW(i).box;
        a=inf(Box(1));
        b=sup(Box(1));
        c=inf(Box(2));
        d=sup(Box(2));
        if a==b && c==d %box is only a point
            plot(a,c, 'LineStyle','none','Marker','o', 'MarkerEdgeColor', ListeW_Farbe);
        elseif a==b || c==d %box is 1-dimensional
            %plot([a,a,b,b,a], [c,d,d,c,c], 'LineWidth',3, 'Marker','o','Color', ListeW_Farbe);
            if a==b
                fill([a-0.1,a-0.1,b+0.1,b+0.1,a-0.1], [c,d,d,c,c], ListeW_Farbe);
            else
                fill([a,a,b,b,a], [c-0.1,d+0.1,d+0.1,c-0.1,c-0.1], ListeW_Farbe);
            end
        else
            fill([a,a,b,b,a], [c,d,d,c,c], ListeW_Farbe);
        end
    end
    
    for i = 1:length(ListNS)
        Box=ListNS(i).box;
        a=inf(Box(1));
        b=sup(Box(1));
        c=inf(Box(2));
        d=sup(Box(2));
        if a==b && c==d %box is only a point
            plot(a,c, 'LineStyle','none','Marker','o', 'MarkerEdgeColor', ListNS(i).Farbe);
        elseif a==b || c==d %box is 1-dimensional
            %plot([a,a,b,b,a], [c,d,d,c,c], 'LineWidth',3, 'Marker','o','Color', ListNS(i).Farbe);
            if a==b
                fill([a-0.1,a-0.1,b+0.1,b+0.1,a-0.1], [c,d,d,c,c], ListNS(i).Farbe);
            else
                fill([a,a,b,b,a], [c-0.1,d+0.1,d+0.1,c-0.1,c-0.1], ListNS(i).Farbe);
            end
        else
            fill([a,a,b,b,a], [c,d,d,c,c], ListNS(i).Farbe);
        end
    end
    
    for i = 1:length(ListS)
        Box=ListS(i).box;
        a=inf(Box(1));
        b=sup(Box(1));
        c=inf(Box(2));
        d=sup(Box(2));
        if a==b && c==d %box is only a point
            plot(a,c, 'LineStyle','none','Marker','o', 'MarkerEdgeColor', ListeS_Farbe);
        elseif a==b || c==d %box is 1-dimensional
            %plot([a,a,b,b,a], [c,d,d,c,c],'LineWidth',3, 'Marker','o', 'Color', ListeS_Farbe);
            if a==b
                fill([a-0.1,a-0.1,b+0.1,b+0.1,a-0.1], [c,d,d,c,c], ListeS_Farbe);
            else
                fill([a,a,b,b,a], [c-0.1,d+0.1,d+0.1,c-0.1,c-0.1], ListeS_Farbe);
            end
        else
            fill([a,a,b,b,a], [c,d,d,c,c], ListeS_Farbe);
        end        
    end
    drawnow;
    hold off;
end
end
