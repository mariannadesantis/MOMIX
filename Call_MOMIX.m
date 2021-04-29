 function ExpiredTimeFlag=Call_MOMIX(optproblem,methodbisect,callMIP,drawCurrentPlots,saveResults,parameter,timel)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOMIX: a branch-and-bound methods for solving Multiobjective Mixed 
%        Integer Convex Problems presented in the paper
%
% [1] "Solving Multiobjective Mixed Integer Convex Optimization Problems"
%     by M. De Santis, G. Eichfelder, J. Niebling, S. Rocktaeschel
%     SIAM Journal on Optimization, 30(4), 3122-3145, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %   Input parameters:
% %   * optproblem: 'T*' (* = 1,...,6)
% %   *methodbisect: {0,1} with br1 = 0; br 2 = 1
% %                 Choose the branching rule (which strategy should be used to bisect a box)
% %                 set methodbisect=0 to choose strategy (br1) detailed in [1]
% %                 set methodbisect=1 to choose strategy (br2) detailed in [1]
% %   *callMIP: {0,1} if set to 0 we call MOMIX_light (NO GUROBI CALL);
% %             Decide if you want enrich the lower bound by computing a further hyperplane
% %             by addressing a single-objective mixed integer problem.
% %             By setting callMIPsolver=1, you run MOMIX 
% %             By setting callMIPsolver=0, you run MOMIX_light
% %   *drawCurrentPlots: {0,1} if set to 0 we do not get any plot along the run
% %   *saveResults: {0,1} if set to 1, results will be saved (and drawCurrentPlots is automatically set to 0)
% %   *parameter: it depends on the test instance - it is used only when
% %               scalable instances are run, check the test instances, to set this parameter properly
% %   *timel: time limit in seconds

%% Initialization
close all;
clear func
clc;
format long;
rng(03677,'twister')
addpath('savefiles')
addpath('test instances')
tic;

global callMIPsolver
callMIPsolver=callMIP;

global timelimit stopAlgo
timelimit=timel;
stopAlgo=0;

if saveResults 
    drawCurrentPlots=0;
end

TI = str2func(optproblem);
if contains(optproblem,'T4')
    m=parameter(1); %number of continuous variables 
    n=parameter(2); %number of integer variables
    [X0,Y0,p,g0,delta,epsilon,epsiloncomp, Q,c,Acon,bcon,Qcon,qcon,ccon]=TI(m,n);
elseif contains(optproblem,'T2') || contains(optproblem,'T3')
    scalablenum=parameter;
    [X0,Y0,p,g0,delta,epsilon,epsiloncomp, Q,c,Acon,bcon,Qcon,qcon,ccon]=TI(scalablenum);
else
    [X0,Y0,p,g0,delta,epsilon,epsiloncomp, Q,c,Acon,bcon,Qcon,qcon,ccon]=TI();
end

if isempty(Q) && isempty(c)
    callMIPsolver=0;
end

%number of continous variables
m=length(X0);
%number of integer variables
n=length(Y0);
%dimension of the preimage space
r=m+n;
%box constraints for the preimage space
Omega0=infsup([inf(X0),inf(Y0)],[sup(X0),sup(Y0)]);
%number of constraints
q=length(g0(mid(Omega0)));

%create folder to save results
if saveResults
    foldername=strcat(char(datetime('now','Format','yyyy_MM_dd_HH_mm')),'_',optproblem,'_branchingrule', num2str(methodbisect),'_MIPsolve', num2str(callMIPsolver), '_numCont', num2str(m), '_numInt', num2str(n));
    dirct=strcat('savefiles/',foldername);
    mkdir(dirct);
    filename = fullfile(dirct, 'results.txt');
    fileID = fopen(filename, 'w');
    fprintf(fileID, ['Objective function: ', optproblem,  ...
                     '\n Numbers of continuous variables           ', num2str(m), ...
                     '\n Number of integer variables               ', num2str(n),...
                     '\n Branching rule                            ', num2str(methodbisect), ...
                     '\n Use MIPsolver to get integer hyperplanes? ', num2str(callMIPsolver), ...
                     '\n delta                                     ', num2str(delta), ...
                     '\n set time limit                            ', num2str(timelimit)]);
end

% objective function variable
global func

% feasibility variable
clear feas
global feas
feas=all(sup(g0(Omega0))<=0);

%initialize lists
ListS_1 = [];
ListNS = [];
%After termination of the algorithm the boxes on
%   -ListS will cover the set of efficient points of the test instance.
%	-ListNS will not contain any efficient point.

%ListPNS = list for points of a potentially nondominated set and their preimages
ListPNS=[];
A=[];

%initialize the list of local upper bounds (ListLUB)
M=max(sup(f(Omega0,func,p)))+1e-6;
Mv=M*ones(1,p);
mv=(min(inf(f(Omega0,func,p)))-1e-6)*ones(1,p);
ListLUB=Mv;


%standardize g0 on Omega0
UBNORMG=norm(sup(abs(g0(Omega0))));
if UBNORMG>0
    g=@(x)g0(x)/UBNORMG;
else
    g=g0;
end

global countLBP GurobiTime fminconTime countRepeatFmincon
countLBP=0;
GurobiTime=0;
fminconTime=0;
countRepeatFmincon=0;

%ListW contains boxes, that still have to be examined, because
%they could not be rejected and do not fulfill the convergence criteria.
ListW = struct('box', Omega0, 'feas',feas,'omegafeas',[],'idealpoint',[]);

%calculating some image points for plots
global FV
FV=[];
 
t_init=toc;

%% MAIN LOOP
tic;
%iteration counter
ind_firstloop=0;

while ~isempty(ListW) && toc<timelimit && ~stopAlgo
    ind_firstloop=ind_firstloop+1;
    %Consider the first box on ListW. This one has the lexicographical
    %smallest ideal point.
    Omega=ListW(1).box;
    %Check if a sufficient condition for the feasibility
    %of all points in this box is fulfilled.
    feas=ListW(1).feas;
    %If a feasible point of MOMINLP in this box is known, use this information.
    omegapotfeas=ListW(1).omegafeas;
    %Delete the current box from ListW.
    ListW=ListW(2:length(ListW));
    %Split the current box into two subboxes.
    if methodbisect==1  % use (br2)
        DIAM=diam(Omega);
        [~,maxdiamidx]=max(DIAM(end:-1:1));
        maxdiamidx=r+1-maxdiamidx;
    else % use (br1)
        DIAM=diam(Omega(m+1:end));
        if sum(DIAM)==0 %all integer variables are fixed
            DIAM=diam(Omega(1:m)); %split at longest continuous edge
            [~,maxdiamidx]=max(DIAM);
        else
            [~,maxdiamidx]=max(DIAM);
            maxdiamidx=m+maxdiamidx;
        end
    end
    Omega1 = Omega;
    Omega2 = Omega;
    if maxdiamidx<=m
        Omega1(maxdiamidx) = infsup(inf(Omega(maxdiamidx)),mid(Omega(maxdiamidx)));
        Omega2(maxdiamidx) = infsup(mid(Omega(maxdiamidx)),sup(Omega(maxdiamidx)));
    elseif mod(inf(Omega(maxdiamidx))+sup(Omega(maxdiamidx)),2)
        Omega1(maxdiamidx) = infsup(ceil(inf(Omega(maxdiamidx))),floor(mid(Omega(maxdiamidx))));
        Omega2(maxdiamidx) = infsup(ceil(mid(Omega(maxdiamidx))),floor(sup(Omega(maxdiamidx))));
    else
        Omega1(maxdiamidx) = infsup(ceil(inf(Omega(maxdiamidx))),floor(mid(Omega(maxdiamidx))));
        Omega2(maxdiamidx) = infsup(ceil(mid(Omega(maxdiamidx)))+1,floor(sup(Omega(maxdiamidx))));
    end  
    %Check if the subboxes can be pruned or fulfill stopping criteria.
    [ListPNS, ListLUB, ListS_1, ListNS, ListW]=Boxcheck2(Omega1,m,n,g0,p, func, ListPNS, ListLUB, ListS_1, ListNS, ListW,epsiloncomp,delta,Q,c,Acon,bcon,Qcon,qcon,ccon);
    [ListPNS, ListLUB, ListS_1, ListNS, ListW]=Boxcheck2(Omega2,m,n,g0,p, func, ListPNS, ListLUB, ListS_1, ListNS, ListW,epsiloncomp,delta,Q,c,Acon,bcon,Qcon,qcon,ccon);
    if drawCurrentPlots
        if r==2
            plot_boxes(Omega0,ListS_1,ListNS,ListW);
        elseif p==2 || p==3
            plot_image(p,FV,ListPNS);
        end    
    end
%     fprintf(['current time for main loop (iteration ', num2str(ind_firstloop), ') is ', num2str(toc), ' s\n']);
end
t1=toc;

if isempty(ListW)
    ExpiredTimeFlag=0;
    h1=plot_boxes(Omega0,ListS_1,ListNS,ListW);

    if saveResults && ~isempty(h1)
        filename=strcat('Boxes');
        saveas(h1,strcat(dirct,'/',filename,'.fig'));
        %saveas(h1,strcat(dirct,'/',filename,'.png'));
    end

    feas=false;
    tic;

    %% Postprocessing
    ListS_2=[];
    ind_secondloop=0;
    while ~isempty(ListS_1)
        ind_secondloop=ind_secondloop+1;
        Omega=ListS_1(1).box;
        H=ListS_1(1).hyper_planes;
        omegafeas=ListS_1(1).omegafeas;
        if isempty(omegafeas)
            omegafeas=mid(Omega);
        end
        a=ListS_1(1).idealpoint;
        ListS_1=ListS_1(2:length(ListS_1));
        if ~isempty(H)
            D=true;
            for s=1:size(ListLUB,1)
                lub=ListLUB(s,:);        
                const=zeros(1,size(H,1));
                for l=1:size(H,1)
                    const(l)=H(l,p+1:end-1)*H(l,1:p)';
                end
                if all(lub*(H(:,1:p)')>= const)  %lub is inside of outer approximation
                    D=false;
                    break;
                end
            end
        else %no hyperplanes are computed, i.e., box is a point
            fomega = f(omegafeas,func,p);
            if any(all((repmat(fomega,size(ListPNS,1),1)>=ListPNS(:,r+1:end)-epsiloncomp)')) && ~any(all((repmat(fomega,size(ListPNS,1),1)==ListPNS(:,r+1:end))'))
                D=true;
            else
                D=false;
            end      
        end    

        if ~D
            ListS_2(end+1).box=Omega;
            ListS_2(end).hyper_planes=H;
            ListS_2(end).omegafeas=omegafeas;
            ListS_2(end).idealpoint=a;
        else
            ListNS(end+1).box=Omega;
            ListNS(end).hyper_planes=H;
            ListNS(end).flag=2;
            ListNS(end).idealpoint=a;

        end
    end
    t2=toc;
    
    %h2=plot_boxes(Omega0,ListS_2,ListNS,ListW);
    %if saveResults && ~isempty(h2)
       % filename=strcat('Boxes2');
       % saveas(h2,strcat(dirct,'/',filename,'.fig'));
       % saveas(h2,strcat(dirct,'/',filename,'.png'));
    %end    

    A=[A;ListPNS(:,1:r)];

    %% Visualization of the results in criterion space
    F=zeros(size(A,1),m);
    for i=1:size(A,1)
        for s=1:p
            F(i,s)=func(A(i,:),s);
        end
    end

    switch p
        case 2
            h4=figure;
            hold on
            if ~isempty(FV)
                plot(FV(:,end-1),FV(:,end),'LineStyle','none','Marker','.','MarkerEdgeColor',0.7*[1 1 1])
            end
            plot(F(:,1),F(:,2),'LineStyle','none','Marker','.','MarkerEdgeColor',[0 0 0])
            xlabel('f1(x)');
            ylabel('f2(x)');
            set(gca,'Color',[0.9 0.9 0.9], 'Fontsize', 15)
        case 3
            if ~isempty(FV)
                h4=scatter3(FV(:,r+1),FV(:,r+2), FV(:,r+3),'.','MarkerEdgeColor', 0.5*[1 1 1]);
                hold on;
                scatter3(F(:,1),F(:,2), F(:,3),'.','MarkerEdgeColor',[0 0 0])
            else 
                h4=scatter3(F(:,1),F(:,2), F(:,3),'.','MarkerEdgeColor',[0 0 0]);
            end        
            xlabel('f1(x)');
            ylabel('f2(x)');
            zlabel('f3(x)');
            set(gca,'Fontsize', 15)    
    end

    if saveResults && ~isempty(h4)
        filename=strcat('Image');
        saveas(h4,strcat(dirct,'/',filename,'.fig'));
        %saveas(h4,strcat(dirct,'/',filename,'.png'));
    end

    clear h1;
    clear h2;
    clear h4;

    if saveResults 
        fprintf(fileID, ['\nResults:',...
                         '\n Numbers of bisections (main loop)                   ', num2str(ind_firstloop), ...
                         '\n Numbers of nodes (main loop)                        ', num2str(2*ind_firstloop+1), ...
                         '\n Computational time for main loop                    ', num2str(t1), 's', ...
                         '\n Numbers of iterations for postprocessing (2nd loop) ', num2str(ind_secondloop), ...
                         '\n Computational time postprocessing                   ', num2str(t2), 's',...
                         '\n Number of found box is L_S                          ', num2str(size(ListS_2,2)),...
                         '\n Computational time for Gurobi                       ', num2str(GurobiTime), ...
                         '\n Computational time for fmincon                      ', num2str(fminconTime), ...
                         '\n Times fmincon was repeated because of bad exitflag  ', num2str(countRepeatFmincon)]);
        fclose(fileID);
        save(strcat(dirct,'/workspace.mat'));
    end
else %computational time expired
    ExpiredTimeFlag=1;
    h1=plot_boxes(Omega0,ListS_1,ListNS,ListW);
    h4=plot_image(p,FV,ListPNS);
    
    if saveResults && ~isempty(h1)
        filename=strcat('Boxes');
        saveas(h1,strcat(dirct,'/',filename,'.fig'));
        %saveas(h1,strcat(dirct,'/',filename,'.png'));
    end
    
    if saveResults &&~isempty(h4)
        filename=strcat('Image');
        saveas(h4,strcat(dirct,'/',filename,'.fig'));
        %saveas(h4,strcat(dirct,'/',filename,'.png'));
    end
    
    clear h1;
    clear h4;
    
    if saveResults
        fprintf(fileID,['\n\nComputational time expired! \n Algorithm could perform ', num2str(ind_firstloop), ' bisections']);
        save(strcat(dirct,'/workspace.mat'));
    end       
        
end
