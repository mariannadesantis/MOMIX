function [ListPNS, ListLUB, ListS_1, ListNS, ListW]=Boxcheck2(Omega,m,n,g0,p,func,ListPNS,ListLUB,ListS_1,ListNS,ListW,epsiloncomp,delta,Q,c,Acon,bcon,Qcon,qcon,ccon)
%% Check if the input boxes can be discarded
%Refer to the paper 
%
% [1] "Solving Multiobjective Mixed Integer Convex Optimization Problems"
%     by M. De Santis, G. Eichfelder, J. Niebling, S. Rocktaeschel
%     SIAM Journal on Optimization, 30(4), 3122-3145, 2020

global callMIPsolver 
global countLBP fminconTime sortTime countRepeatFmincon timelimit stopAlgo
global FV

box_infeas=false;

%standardize g on Omega
UBNORMG=norm(sup(abs(g0(Omega))));
if UBNORMG>0
    g=@(x)g0(x)/UBNORMG;
else
    g=g0;
end
r=m+n;

% feas indicates whether the superbox of Omega fulfill a sufficient
% condition for containing only feasible points of the relaxed problem
global feas
hilf=feas;

infeas=false;
if ~infeas
    countLBP=countLBP+1;
    omegastart=mid(Omega);  
    %% Apply algorithm 2 in [1]
    Omegal=inf(Omega);
    Omegau=sup(Omega);
    %Compute the ideal point of the continuous relaxation of the problem
    if all(Omegal==Omegau) % The box consists of only one point, it happens only if all variables are integer
        omegafeas=Omegal;
        if (feas || all(g(omegafeas)<=0))
            FV=[FV; omegafeas,f(omegafeas,func,p)];
        end
        fomega = f(omegafeas,func,p);
        a=fomega;
        [ListPNS, add_onepoint] = updateListPNS(Omegal,m, func,p,g, ListPNS, epsiloncomp);
        H=[];
        if add_onepoint || any(all((repmat(fomega,size(ListPNS,1),1)==ListPNS(:,r+1:end))'))
            D=false;
            if add_onepoint
                ListLUB=lub3(ListLUB,fomega,p);
            end
        else
            D=true;
        end
    else
        a=zeros(1,p);
        add=zeros(1,p);
        fmincopt=optimset('Algorithm','interior-point','Display','off', 'MaxIter', 10000);
        % Address the single-objective continuous convex problems to
        % compute the ideal point:
        for i=1:p
            func_i=@(omega)func(omega,i);
            fminconstart=toc;
            [omega,a(i),exitflag] = fmincon(func_i,omegastart,[],[],[],[],Omegal,Omegau,@(omega)constr(g,omega,m),fmincopt);
            fminconTime=fminconTime+toc-fminconstart;
            while toc<timelimit && exitflag~=1 && exitflag~=2 && exitflag~=-2 
                countRepeatFmincon=countRepeatFmincon+1;
                fminconstart=toc;
                [omega,a(i),exitflag] = fmincon(func_i,omega,[],[],[],[],Omegal,Omegau,@(omega)constr(g,omega,m),fmincopt);
                fminconTime=fminconTime+toc-fminconstart;
            end
            if toc>= timelimit
                stopAlgo=1;
                return
            end
            %If fmincon terminates with exitflag=-2, no feasible points for 
            % the continuous relaxation were found. Hence, Omega does not 
            % contain any feasible points for MOMINLP and can be discarded 
            % without further investigation.
            if exitflag==-2
                box_infeas=true;
                %fprintf('convex problem is infeasible\n');
                break;
            else
                omegaround=[omega(1:m),round(omega(m+1:end))];
                if (feas || all(g(omegaround)<=0))
                    FV=[FV; omegaround, f(omegaround,func,p)];
                end
                omegafeas=omegaround;
                % Add rounded omega and its func value to ListPNS, if it is
                % feasible for MOMINLP
                [ListPNS, add(i)] = updateListPNS(omegaround,m, func,p,g, ListPNS, epsiloncomp);
                % If omegaround and its func value are added to ListPNS, we
                % update the LUB.
                if add(i)
                    ListLUB=lub3(ListLUB,f(omegaround,func,p),p);                    
                end
            end
        end
   
        if any(add)
            box_infeas=false;
        end
        if ~box_infeas
            % Initialize supporting hyperplanes. last column: indicator
            % which type of hyperplane: init: 0; convex: 1; integer: 2
            H=[eye(p),repmat(a,p,1),zeros(p,1)];
            % indicator for discarding Omega
            D=true;
            % objective function for P_lub
            fun=@(omega)omega(r+1);
            ListLUBfix=ListLUB;

            for s=1:size(ListLUBfix,1)
                lub=ListLUBfix(s,:);
                const=zeros(1,size(H,1));
                for j=1:size(H,1)
                    const(j)=H(j,p+1:end-1)*H(j,1:p)';
                end

                if ~isempty(H) && all(lub*(H(:,1:p)')>= const) % if lub is inside the outer approximation
                    % starting point for P_lub
                    omega0=[omegastart,max(f(omegastart,func,p)-lub)+1];
                    % Solve P_lub.
                    fminconstart=toc;
                    [z,~,exitflag,~,lagrange] = fmincon(fun,omega0,[],[],[],[],[inf(Omega),-Inf],[sup(Omega),Inf],@(z)P_lubconstr(func,z,g,lub,p),fmincopt);
                    fminconTime=fminconTime+toc-fminconstart;
                    %fprintf(['exitflag for a OP to compute (relaxed) hyperplane is ', num2str(exitflag), '\n']);
                    while toc<timelimit && exitflag~=1 && exitflag~=2 && exitflag~=-2
                        %fmincopt=optimset('Algorithm','interior-point','Display','off', 'MaxIter', 10000);
                        countRepeatFmincon=countRepeatFmincon+1;
                        fminconstart=toc;
                        [z,~,exitflag,~,lagrange] = fmincon(fun,z,[],[],[],[],[inf(Omega),-Inf],[sup(Omega),Inf],@(z)P_lubconstr(func,z,g,lub,p),fmincopt);
                        fminconTime=fminconTime+toc-fminconstart;
                        %fprintf(['exitflag for a repeated convex OP is ', num2str(exitflag), '\n']);
                    end
                    if toc>= timelimit
                        stopAlgo=1;
                        return
                    end
                    omega_opt=z(1:r);
                    t_opt=z(r+1);
                    % Add rounded omega_opt and its func value to ListPNS, 
                    % if it is feasible for MOMINLP
                    omegaround=[omega_opt(1:m),round(omega_opt(m+1:end))];
                    omegafeas=omegaround;
                    if (feas || all(g(omegaround)<=0))
                        FV=[FV;omegaround, f(omegaround,func,p)];
                    end
                    [ListPNS, add(s)] = updateListPNS(omegaround,m, func,p,g, ListPNS, epsiloncomp);
                    % If omegaround and its func value are added to ListPNS, 
                    % we update the LUB.
                    if add(s)
                        ListLUB=lub3(ListLUB,f(omegaround,func,p),p);                        
                    end

                    %normal vector of the hyperplane regarding lub
                    lambda=lagrange.ineqnonlin(1:p);
                    hp= lub+t_opt*ones(1,p);
                    if t_opt<=epsiloncomp && callMIPsolver % Lub lies above the convex hyperplane
                        [hp_int,ListPNS,ListLUB,exitflag,omegafeas]=integerHyperplanes(lambda,Omega,func,p,m,n,g, ListPNS, ListLUB, epsiloncomp,Q,c,Acon,bcon,Qcon,qcon,ccon);
                        %fprintf(['t*= ', num2str(t_opt), ' -> integer problem is solved with exitflag ', num2str(exitflag), '\n']);
                        if exitflag==-2
                            box_infeas=true;
                            %fprintf('integer problem is infeasible\n');
                            break;
                        end
                        if (feas || all(g(omegafeas)<=0))
                            FV=[FV; omegafeas, f(omegafeas,func,p)];
                        end
                        H=[H;lambda',hp_int,2];
                        if lub*lambda>=hp_int*lambda 
                            % lub lies above integer hyperplane and then 
                            % the box should be bisected
                            D=false;
                            break;
                        end
                    elseif t_opt<=epsiloncomp %lub lies above convex hyperplane and MIPsolver should not be called
                        H=[H;lambda',hp,1];
                        D=false;
                    else % topt is positive
                        H=[H;lambda',hp,1];
                        % fprintf('lub lies under convex hyperplane\n');
                    end
                end
            end
        end
    end
else
    % in this case, Omega does not contain any feasible points for 
    % the continuous relaxation
    % Hence, it does not contain any feasible point for MOMINLP.
    box_infeas=true;
    %fprintf('intervalarithmetic on g is greater than 0\n');
end
%% Update lists.

if box_infeas
    ListNS(end+1).box=Omega;
    ListNS(end).flag=0;
%     fprintf('Box discarded\n');
else
    if D
        ListNS(end+1).box=Omega;
        ListNS(end).hyper_planes=H;
        ListNS(end).idealpoint=a;
        ListNS(end).flag=1;
%         fprintf('Box discarded\n');
    elseif norm(sup(Omega)-inf(Omega))<delta %Termination Rule is satisfied
        ListS_1(end+1).box=Omega;
        ListS_1(end).hyper_planes=H;
        if all(round(omegafeas(m+1:end))==omegafeas(m+1:end)) &&( feas || all(g(omegafeas)<=0) )
            ListS_1(end).omegafeas=omegafeas;
        else
            ListS_1(end).omegafeas=[];
        end

        ListS_1(end).idealpoint=a;
%         fprintf('Box in solution list\n');
    else
        ListW(end+1).box=Omega;
        ListW(end).feas=feas;
        if all(round(omegafeas(m+1:end))==omegafeas(m+1:end)) &&( feas || all(g(omegafeas)<=0) )
            ListW(end).omegafeas=omegafeas;
        else
            ListW(end).omegafeas=[];
        end
        ListW(end).idealpoint=a;
%         fprintf('Box in working list\n');
        if length(ListW)>1
            sortbegin=toc;
            ListW=sort_list(ListW);
            sortTime=sortTime+toc-sortbegin;
        end
    end
end

%Reset feas for next box.
feas=hilf;
end