function minlp = PaDOA_subalg_setup(solver_struct,opts,discrete)
    import casadi.*
    opts.MINLPsolver_opts.discrete = discrete;
    if strcmp(opts.MINLP_solver,'bonmin')
        minlp = nlpsol('solver','bonmin',solver_struct,opts.MINLPsolver_opts);
    elseif strcmp(opts.MINLP_solver,'gurobi')
        minlp = qpsol('solver','gurobi',solver_struct,opts.MINLPsolver_opts);
    elseif strcmp(opts.MINLP_solver,'cplex')
        minlp = qpsol('solver','cplex',solver_struct,opts.MINLPsolver_opts);
    elseif strcmp(opts.MINLP_solver,'ipopt')
        if sum(discrete)~=0
           error('IPOPT can only be applied to real-valued problems. Non-zero discrete vector detected.') 
        end
        minlp = nlpsol('solver','ipopt',solver_struct,opts.MINLPsolver_opts);
    end
end