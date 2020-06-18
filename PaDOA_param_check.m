function [update_hyper_planes_with_U,remove_inac_hyperplanes,limited_hyperplane_gen,upper_bound_MILP_epigraph,lower_bound_MILP_epigraph,max_hyperplanes,Us,Ls,eps,maxiter] = PaDOA_param_check(params,N)
    if isfield(params,'update_hyper_planes_with_U')
        update_hyper_planes_with_U = params.update_hyper_planes_with_U;
    else
        update_hyper_planes_with_U = 1;
    end
    if isfield(params,'remove_inac_hyperplanes')
        remove_inac_hyperplanes = params.remove_inac_hyperplanes;
    else remove_inac_hyperplanes = 1;
    end
    if isfield(params,'limited_hyperplane_gen')
        limited_hyperplane_gen = params.limited_hyperplane_gen;
    else limited_hyperplane_gen = 0;
    end
    if isfield(params,'upper_bound_MILP_epigraph')
        upper_bound_MILP_epigraph = params.upper_bound_MILP_epigraph;
    else upper_bound_MILP_epigraph = 0;
    end
    if isfield(params,'lower_bound_MILP_epigraph')
        lower_bound_MILP_epigraph = params.lower_bound_MILP_epigraph;
    else lower_bound_MILP_epigraph = 1;
    end
    if isfield(params,'max_hyperplanes')
        max_hyperplanes = params.max_hyperplanes;
    else max_hyperplanes = inf;
    end
    if isfield(params,'Us')
        Us = params.Us;
    else Us = inf(N,1);
    end
    if isfield(params,'Ls')
        Ls = params.Ls;
    else Ls = -inf(N,1);
    end
    if isfield(params,'eps')
        eps = params.eps;
    else eps = 10^-3;
    end
    if isfield(params,'maxiter')
        maxiter = params.maxiter;
    else maxiter = 100;
    end
end
