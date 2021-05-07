# MOMIX
a branch and bound method for multi-objective mixed integer convex optimization problems

The present repository contains the implementation of MOMIX and MOMIX_light:
the branch-and-bound methods for solving Multiobjective Mixed Integer Convex
Problems presented in the paper

[1] "Solving Multiobjective Mixed Integer Convex Optimization Problems"
by M. De Santis, G. Eichfelder, J. Niebling, S. Rocktaeschel
SIAM Journal on Optimization, 30(4), 3122-3145, 2020.

In case you want to reproduce the results on your own machine,
please install GUROBI and INTLAB on your PC


=========================================================================================

The repository provides

1) the test instances (as .m files) used for the numerical 
   experiments in [1]:(T1),(T2),(T3),(T4),(T5),(T6).

2) the following .m files for running MOMIX:

- start_session.m : in order to start the code, it is necessary to start INTLAB and GUROBI
                    change the directories according to give the proper links

- Testrun.m : example of call for MOMIX, to launch the numerical experiments of paper [1];

- Call_MOMIX.m : main file
  (take care of the folder directories according to your operation system)
  Call_MOMIX.m has the following input parameters:

  * optproblem : 'T*' (* = 1,...,6)
  
  * methodbisect : {0,1} choose the branching strategy: for (br1) choose 0; for (br2) choose 1

  * callMIP : {0,1} if set to 0 we call MOMIX_light (NO GUROBI CALL, i.e. No integer hyperplanes);

  * drawCurrentPlots : {0,1} if set to 0 we do not get any plot along the run of the algorithm

  * saveResults : {0,1} if set to 1, results will be saved (and drawCurrentPlots is automatically set to 0)

  * parameter : it depends on the test instance - it is used only when scalable instances are run
              - set according to what is mentioned within the instance files 

  * timel : time limit in seconds


- Boxcheck2.m : Discarding procedure (Algorithm 2 in [1])

- sort_list.m : it sorts the boxes within the working list, in order to choose the box with lexicographic smallest ideal point

- updateListPNS.m : it updates the set L_PNS

- lub3.m : procedure to compute the local upper bound set from the L_PNS list

- finde.m : called by lub3.m

- integerHyperplanes.m : called by Boxcheck2.m in case an integer hyperplane should be computed

- GurobiCall.m : called by integerHyperplanes.m  to solve the MIQP subproblems in the bounding procedure

- f.m : to evaluate the objective functions

MATLAB files needed to call FMINCON:
- constr.m : to evaluate the constraints (input for fmincon)
- P_lubconstr.m : needed to call fmincon

MATLAB files for the plot:
- plot_boxes.m : to plot boxes for 2D instances
- plot_image.m : to plot L_PNS and the images of all feasible points found along the run of the algorithm (2D and 3D)


...enjoy, buon divertimento, viel Spass :-)

=============================================================================================
