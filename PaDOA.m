function [xopt, logg, tot_time, PADOA_obj] = PaDOA( obj_funs,cons_funs,lb,ub,A,B,c,x0,z0,params,opts)

import casadi.*

good_statuses = {'integer optimal, tolerance','OPTIMAL','SUCCESS','integer optimal solution'};                                
milp_opts = opts.milp_opts;

% number of subproblems
N  = length(obj_funs);

%% built local subproblems and CasADi functions
all_real_lb = [];
all_real_ub = [];
remove_x = zeros(1,N);
for i=1:N
    if length([x0{i};z0{i}])==1 %CasADi has problems with single variable problems :(
        x_cas = SX.sym('temp',2,1);
        obj_funs{i} = Function('f',{[x_cas]},{obj_funs{i}(x_cas(1))});
        cons_funs{i} = Function('f',{[x_cas]},{cons_funs{i}(x_cas(1))});
        lb{i} = [0;lb{i}];
        ub{i} = [0;ub{i}];
        x0{i} = [0;x0{i}];
        A{i}  = [0,A{i}];
        remove_x(i)=1;
    end
    
    % real-valued vars
    rvars(i)    = length(x0{i});
    xiCas{i}   = SX.sym(strcat('xi_',num2str(i)),rvars(i),1); 
    % local integers (parameterized in MINLP)
    ivars(i)    = length(z0{i});
    zCas{i}    = SX.sym(strcat('zi_',num2str(i)),ivars(i),1); 
    %re-index bounds
    all_real_lb = [all_real_lb;lb{i}(1:length(x0{i}))];
    all_real_ub = [all_real_ub;ub{i}(1:length(x0{i}))];
end

for i=1:N
    % integer variables
    zetaCas = SX.sym('zeta',ivars(i),1); 
    discrete{i} = [zeros(1,rvars(i)),ones(1,ivars(i))];
    discrete_minlp{i} = [zeros(1,sum(rvars)),ones(1,ivars(i))];
    
    % add bounds on all real vars and local int vars
    lb_minlp{i} = [all_real_lb;lb{i}(length(x0{i})+1:length(x0{i})+length(z0{i}))];
    ub_minlp{i} = [all_real_ub;ub{i}(length(x0{i})+1:length(x0{i})+length(z0{i}))];
    
    % non-local ints parameterized
    pCas        = [vertcat(zCas{1:i-1});vertcat(zCas{i+1:N})];
                
    % objective function for local MINLP's
    obj_Loc_Cas = obj_funs{i}([xiCas{i};zetaCas]);
    index = 1:N; index(i)=[];
    for j = index
        obj_Loc_Cas = obj_Loc_Cas + obj_funs{j}([xiCas{j};zCas{j}]);
    end
    
    % local inequality constraints
    cons_Cas  = cons_funs{i}([xiCas{i};zetaCas]);
    for j = index
        cons_Cas = [cons_Cas; cons_funs{j}([xiCas{j};zCas{j}])];
    end
    cons_Cas = [cons_Cas;horzcat(A{:})*vertcat(xiCas{:}) + horzcat(B{:})*[vertcat(zCas{1:i-1}); zetaCas; vertcat(zCas{i+1:N})]-c];
    cons_Cas = [cons_Cas;-horzcat(A{:})*vertcat(xiCas{:}) - horzcat(B{:})*[vertcat(zCas{1:i-1}); zetaCas; vertcat(zCas{i+1:N})]+c];
    nCons{i} = length(cons_Cas);

    %Gradient of local objective wrt x
    alphaCas    = gradient(obj_funs{i}([xiCas{i};zetaCas]),xiCas{i});
    alpha_fun{i}    = Function(['alpha' num2str(i)],{xiCas{i},zetaCas},{alphaCas});
    
    %Gradient of local objective wrt z
    betaCas    = gradient(obj_funs{i}([xiCas{i};zetaCas]),zetaCas);
    beta_fun{i}    = Function(['beta' num2str(i)],{xiCas{i},zetaCas},{betaCas});

    % set up local solvers
    milp_opts.discrete = [zeros(1,sum(rvars)),ones(1,sum(ivars)),zeros(1,N)];
    
    minlp_struct = struct('x',[vertcat(xiCas{:});zetaCas],'f',obj_Loc_Cas,'g',cons_Cas,'p',pCas); 
    minlp{i} = PaDOA_subalg_setup(minlp_struct,opts,discrete_minlp{i});        
end

% set up params
[update_hyper_planes_with_U,remove_inac_hyperplanes,limited_hyperplane_gen,upper_bound_MILP_epigraph,lower_bound_MILP_epigraph,max_hyperplanes,Us,Ls,eps,maxiter] = PaDOA_param_check(params,N);

%MILP setup
ell = SX.sym('l',[N,1]);
x_glob = SX.sym('x',[sum(rvars),1]);
z_glob = SX.sym('z',[sum(ivars),1]);
lbounds = vertcat(lb{:});
ubounds = vertcat(ub{:});
lbx = lbounds(find(1-horzcat(discrete{:})));
ubx = ubounds(find(1-horzcat(discrete{:})));
lbz = lbounds(find(horzcat(discrete{:})));
ubz = ubounds(find(horzcat(discrete{:})));

neqs = length(c); % number of consensus constraints
milp_cons = [horzcat(A{:}) horzcat(B{:})]*[x_glob;z_glob]; % add the consensus constraints to the MILP
obj_cent = 0;
for i = 1:N
    rn_ind = 1+sum(rvars(1:i-1)):sum(rvars(1:i-1))+rvars(i);
    in_ind = 1+sum(ivars(1:i-1)):sum(ivars(1:i-1))+ivars(i);
    milp_cons = [milp_cons; cons_funs{i}([x_glob(rn_ind);z_glob(in_ind)])]; % add the local inequality constraints to the MILP
    obj_cent = obj_cent + obj_funs{i}([xiCas{i};zCas{i}]); % construct the overall objective function
end
obj_fun_cent = Function('f',{[vertcat(xiCas{:});vertcat(zCas{:})]},{obj_cent});

%initialization
i   = 1;
U = inf;
if upper_bound_MILP_epigraph
    U = sum(full(Us));
end
tot_time = 0;
for j = 1:N
z{j} = [vertcat(z0{1:j-1});vertcat(z0{j+1:N})];
loc_sol{j} = [zeros(sum(rvars),1);z0{j}];
end
epi_cons_A = sum(ell);
epi_cons_rhs = 9999;
epi_cons_lhs = -9999;
hyperplanes = [];
logg.x_plus = [];
logg.z_plus = [];
logg.stats = [];
logg.time.subtime = [];
logg.time.milptime = [];
logg.time.hyptime = [];
logg.U = [];
logg.L = [];

%begin iterating
while i <= maxiter 
   
   sub_time = zeros(1,N);
   infeas_MINLPS = 0;
    for j=1:N
        obj_val = 0;
tic             
        % solve local MINLP's 
        [loc_sol{j},minlp_status] = PaDOA_MINLP(minlp{j},loc_sol{j},z{j},lb_minlp{j},ub_minlp{j},-inf(nCons{j},1),zeros(nCons{j},1));
        
        % logging
        logg.stats = [logg.stats;minlp_status];
        if ~any(strcmp(good_statuses,minlp_status))
            infeas_MINLPS = infeas_MINLPS + 1;
            disp('Infeasible MINLP detected')
        else       
            %Store local solutions
            rn_ind = 1+sum(rvars(1:j-1)):sum(rvars(1:j-1))+rvars(j);
            in_ind = 1+sum(ivars(1:j-1)):sum(ivars(1:j-1))+ivars(j);

            x_i_opt{j} = loc_sol{j}(rn_ind);
            z_opt{j} = loc_sol{j}(sum(rvars)+1:end);

            % Compute Upper Bound
            real_loc_sol = loc_sol{j}(1:sum(rvars));
            int_loc_sol = [z{j}(1:sum(ivars(1:j-1)));z_opt{j};z{j}(sum(ivars(1:j-1))+1:end)];
            obj_val = obj_fun_cent([real_loc_sol;int_loc_sol]);
            if or(update_hyper_planes_with_U && full(obj_val)<U, ~update_hyper_planes_with_U)
                %Compute Hyperplanes
                if limited_hyperplane_gen
                   hyp_ind=j;
                else
                   hyp_ind=1:N; 
                end
                for n = hyp_ind 
                    loc_x_index = 1+sum(rvars(1:n-1)):sum(rvars(1:n-1))+rvars(n);
                    loc_z_index = 1+sum(ivars(1:n-1)):sum(ivars(1:n-1))+ivars(n);
                    xsol_loc = real_loc_sol(loc_x_index);
                    zsol_loc = int_loc_sol(loc_z_index);

                    [alpha,beta,gamma] = PaDOA_add_HP(xsol_loc,zsol_loc,obj_funs{j},alpha_fun{j},beta_fun{j}); 

                    xn = x_glob(loc_x_index);  
                    zn = z_glob(loc_z_index);  

                    hyperplanes = [hyperplanes; alpha'*xn + beta'*zn + gamma - ell(n)];  
                end
            end
            
             if full(obj_val)<U
                U = full(obj_val);
                if upper_bound_MILP_epigraph
                    % update epigraph upper bound
                    epi_cons_rhs = U;
                end
             end 
        end
sub_time(j) = toc;
    end
tot_time = tot_time + max(sub_time);
logg.time.subtime = [logg.time.subtime,max(sub_time)];

if infeas_MINLPS == N
   warning('All subproblems failed to return a feasible solution.') 
   keyboard
end
if length(hyperplanes)>max_hyperplanes
    hyperplanes(1:length(hyperplanes)-max_hyperplanes)=[];
end
if or(i == 1,remove_inac_hyperplanes==0)
   active_hyp = 1:length(hyperplanes); 
end

% add hyperplanes to MILP
milp_cons_temp = [milp_cons;hyperplanes(active_hyp)];
nineqs = size(milp_cons_temp,1);
nepicons = 0;
if upper_bound_MILP_epigraph
    nepicons = nepicons + 1;
    milp_cons_temp = [milp_cons_temp;epi_cons_A-epi_cons_rhs];
end
if lower_bound_MILP_epigraph
    nepicons = nepicons + 1;
    milp_cons_temp = [milp_cons_temp;-epi_cons_A+epi_cons_lhs];
end

disp(['Iteration: ',num2str(i),', Hyperplanes: ',num2str(length(active_hyp)),', Epigraph constraints: ',num2str(nepicons)])
   
%Solve MILP
tic
    milp =   struct('x',[x_glob;z_glob;ell],'f',sum(ell),'g',milp_cons_temp);
    cas =   qpsol('solver','cplex',milp,milp_opts); 
    [milp_sol,milp_status] = PaDOA_MINLP(cas,[],[],[lbx;lbz;Ls],[ubx;ubz;Us],[c;-inf(nineqs+nepicons-neqs,1)],[c;zeros(nineqs+nepicons-neqs,1)]);

milp_time = toc; 
tot_time = tot_time + milp_time;
logg.time.milptime = [logg.time.milptime,milp_time];

    %check if anything is breaking...
    logg.stats = [logg.stats;string(milp_status)];
    if ~any(strcmp(good_statuses,minlp_status))
        warning('MILP failed. Resolving without epigraph constraints.')
        milp =   struct('x',[x_glob;z_glob;ell],'f',sum(ell),'g',milp_cons(1:nineqs+neqs)); % remove epigraph constraints
        cas =   qpsol('solver','cplex',milp,milp_opts);
        % try solving without epigraph constraints
        [milp_sol,~] = PaDOA_MINLP(cas,[],[],[lbx;lbz;-inf(N,1)],[ubx;ubz;inf(N,1)],[c;-inf(nineqs-neqs,1)],[c;zeros(nineqs-neqs,1)]);
    end
    
    x_plus = milp_sol(1:sum(rvars));
    z_plus = milp_sol(sum(rvars)+1:sum(rvars)+sum(ivars));
    y_plus = milp_sol(sum(rvars)+sum(ivars)+1:end);
    
    if full(sum(y_plus)) < U
       U = full(sum(y_plus)); %update upper bound
    end
    
    logg.U = [logg.U, U];
    logg.L = [logg.L, sum(y_plus)];
    
    if remove_inac_hyperplanes
        % remove inactive hyperplanes
        hyperplane_fun = Function('h',{[x_glob;z_glob;ell]},{hyperplanes});
        h_vals = full(hyperplane_fun(milp_sol));
        n_inac = sum(h_vals<-10^-4);
        if n_inac>0
            active_hyp = find(h_vals>=-10^-4);
        end
    end
    
    if lower_bound_MILP_epigraph
        % update epigraph lower bound
        epi_cons_lhs = sum(full(y_plus));
    end

    for j = 1:N
tic
        %Compute hyperplanes at MILP solution
        loc_x_index = 1+sum(rvars(1:j-1)):sum(rvars(1:j-1))+rvars(j);
        loc_z_index = 1+sum(ivars(1:j-1)):sum(ivars(1:j-1))+ivars(j);
        xsol_loc = x_plus(loc_x_index);
        zsol_loc = z_plus(loc_z_index);
        
        [alpha,beta,gamma] = PaDOA_add_HP(xsol_loc,zsol_loc,obj_funs{j},alpha_fun{j},beta_fun{j});
        
        xn = x_glob(loc_x_index);
        zn = z_glob(loc_z_index);
        hyperplanes = [hyperplanes; alpha'*xn + beta'*zn + gamma - ell(j)]; 
hyp_time(j) = toc;
        %Update integer variables
        z{j} = full(z_plus);
        z{j}(loc_z_index)=[];
    end
    
tot_time = tot_time + max(hyp_time); 
logg.time.hyptime = [logg.time.hyptime,max(hyp_time)];

    %Termination check
    if U - full(sum(y_plus)) < eps
        logg.x_plus  = [logg.x_plus,x_plus];
        logg.z_plus  = [logg.z_plus,z_plus];
        break
    end
    
    %next iteration and track progress
    i = i+1;
    logg.x_plus  = [logg.x_plus,x_plus];
    logg.z_plus  = [logg.z_plus,z_plus];
    
    disp(['Current upper bound: ',num2str(U)])
    disp(['Current lower bound: ',num2str(sum(full(y_plus)))])
    disp(' ')

end

for i = 1:N
    % if a dummy variable was added, then remove it here
   if remove_x(i)==1
       loc_x_index = 1+sum(rvars(1:j-1)):sum(rvars(1:j-1))+rvars(j);
       if length(x_plus(loc_x_index))<2
           x_plus = [];
       else
           x_plus(loc_x_index(1)) = [];
       end
   end
end

% Output results
disp(['Current upper bound: ',num2str(U)])
disp(['Current lower bound: ',num2str(sum(full(y_plus)))])

xopt = [x_plus;z_plus];
PADOA_obj = 0;
for j = 1:N
    loc_x_index = 1+sum(rvars(1:j-1)):sum(rvars(1:j-1))+rvars(j);
    loc_z_index = 1+sum(ivars(1:j-1)):sum(ivars(1:j-1))+ivars(j);
    PADOA_obj = PADOA_obj + full(obj_funs{j}([x_plus(loc_x_index);z_plus(loc_z_index)]));
end
end