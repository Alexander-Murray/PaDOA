# PaDOA
Algorithm for partially distributed outer approximation of MICPs. This code uses the open-source tool CasADi:  
https://web.casadi.org/  

Input problem must be partitioned:
**obj_funs:** cell array of CasADi functions
**cons_funs:** cell array of CasADi functions
**lb:** cell array of vectors  
**ub:** cell array of vectors  
**A:** cell array of matrices  
**B:** cell array of matrices  
**c:** vector  
**x0:** cell array of vectors  
**z0:** cell array of vectors  
**params:** struct
**opts:** struct  

Params will assume default values if left empty. Possible params to set:  
**update_hyper_planes_with_U:** 1,0. Whether or not to add hyperplanes only when an MINLP solution improves the upper bound. Default: 1  
**remove_inac_hyperplanes:** 1,0. Whether or not to check active hyperplanes at the MILP solutio and remove those that are inactive. Default: 1  
**limited_hyperplane_gen:** 1,0. Whether to add 1 instead of N hyperplanes from each MINLP solution. Default: 0  
**upper_bound_MILP_epigraph:** 1,0. Whether or or not to add an upper bound on the MILP. Default: 0  
**lower_bound_MILP_epigraph:** 1,0. Whether or or not to add a lower bound on the MILP. Default: 1  
**max_hyperplanes:** double. Start discarding old hyperplanes after this amount has been added to the MILP. Default: inf  
**Us:** vector. User provided upper bounds on each MINLP subproblem. Default: inf(N,1);   
**Ls:** vector. User provided lower bounds on each MINLP subproblem. Default: -inf(N,1);  
**eps:** double. Termination tolerance. Default: 10^-3;  
**maxiter:** int. Iteration limit. Default: 100;  

Opts must be set. Algorithm options are:  
**MILP_solver:** 'ipopt', 'bonmin', 'gurobi', 'cplex'  
**MINLP_solver:** 'ipopt', 'bonmin', 'gurobi', 'cplex'  
**MINLPsolver_opts:** used to pass options to the MINLP solver   
**milp_opts:** used to pass options to the MILP solver  

The syntax for using PaDOA is shown in the EXAMPLE_MIP.  

To cite use:  
@Article{Murray_2019,  
  author  = {Murray, A. and Faulwasser, T. and Hagenmeyer, V. and Villanueva, M.E. and Houska, B.},  
  title   = {Partially Distributed Outer Approximation},  
  journal = {Journal of Global Optimization},  
  year    = {2020},  
}