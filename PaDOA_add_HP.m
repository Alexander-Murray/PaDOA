function [alpha,beta,gamma] = PaDOA_add_HP(xsol_loc,zsol_loc,obj_funs,alpha_fun,beta_fun)
    %Compute hyperplanes at MILP solution
    alpha = alpha_fun(xsol_loc,zsol_loc);
    beta  = beta_fun(xsol_loc,zsol_loc);
    if isempty(xsol_loc)
        gamma = obj_funs(zsol_loc)-beta.'*zsol_loc;
    elseif isempty(zsol_loc)
        gamma = obj_funs(xsol_loc)-alpha.'*xsol_loc;
    else
        gamma = obj_funs([xsol_loc;zsol_loc])-alpha.'*xsol_loc-beta.'*zsol_loc;
    end
end