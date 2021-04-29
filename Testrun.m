%%% Example for the call of MOMIX:
%%% according to the experiments performed in the paper:
%%%
%%% "Solving Multiobjective Mixed Integer Convex Optimization Problems"
%%% by M. De Santis, G. Eichfelder, J. Niebling, S. Rocktaeschel
%%% SIAM Journal on Optimization, 30(4), 3122-3145, 2020.

clear all;

%start_session;

% Name the test problems
problem{1}='T1';
problem{2}='T2';
problem{3}='T3';
problem{4}='T4';
problem{5}='T5';
problem{6}='T6';

%Call_MOMIX(optproblem,methodbisect,callMIP,drawCurrentPlots,saveResults,parameter,timel)

% First experiment:
Call_MOMIX('T1',0, 0, 0, 1, 0, 1800)

%Second experiment (first line, first column)
%Call_MOMIX('T2',0, 1, 0, 1, 3, 1800)