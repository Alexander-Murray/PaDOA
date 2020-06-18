close all;
clear all;
clc;

import casadi.*

%% Problem parameters
N = 10; % number of partitions
NC = 5; % number of consensus constraints
nx = 10*ones(N,1); % number of real vars in each partition
nz = 10*ones(N,1); % number of discrete vars in each partition
nxz = nx+nz;

%% Problem setup
for i = 1:N
    x = SX.sym(strcat('x',num2str(i),'_'),nx(i),1); % local real vars
    z = SX.sym(strcat('z',num2str(i),'_'),nz(i),1); % local discrete vars
    
    % objective function
    obj{i} = rand(1,nx(i))*x + rand(1,nz(i))*z;
    obj_funs{i} = Function('f',{[x;z]},{obj{i}});
    
    % local constraints
    cons{i} = x.'*x + i*z - 5*rand;
    cons_funs{i} = Function('f',{[x;z]},{cons{i}});
    
    % consensus constraint matrices
    A{i} = rand(NC,nx(i));
    B{i} = zeros(NC,nz(i));
    AB{i} = [A{i}, B{i}];
    
    % initial guess
    x0{i} = zeros(nx(i),1);
    z0{i} = zeros(nz(i),1);
    xz0{i} = zeros(nx(i)+nz(i),1);
    
    % box constraints
    lb{i} = -2*ones(nxz(i),1);
    ub{i} = 2*ones(nxz(i),1);
    
    % vector denoting which variables are integer valued
    discrete{i}=[zeros(1,nx(i)),ones(1,nz(i))];
end
c = zeros(NC,1); % RHS of consensus constraints
lam0 = zeros(NC,1); % initial guess of Lagrange multiplier of consensus constraints

%% Centralized solution

% Centralized problem formulation
x_glob =  SX.sym('x',[sum(nxz),1]);

obj_cent_temp = 0;
cons_cent_temp= horzcat(AB{:})*x_glob-c;
for i = 1:N
    x_part = x_glob(sum(nxz(1:i-1))+1:sum(nxz(1:i)));

    obj_cent_temp = obj_cent_temp + obj_funs{i}(x_part);
    
    cons_cent_temp = [cons_cent_temp;cons_funs{i}(x_part)];
end
obj_cent = Function('f',{x_glob},{obj_cent_temp});

cons_cent = Function('g',{x_glob},{cons_cent_temp});

discrete_cent = horzcat(discrete{:});

lbg = [zeros(NC,1);-Inf(length(cons_cent_temp)-NC,1)];
ubg = [zeros(NC,1);zeros(length(cons_cent_temp)-NC,1)];

% solve with bonmin
tic
cent_opts.verbose_init = 0;
cent_opts.bonmin.print_level = 0;
cent_opts.bonmin.bb_log_level = 0;
cent_opts.bonmin.fp_log_level = 0;
cent_opts.bonmin.lp_log_level = 0;
cent_opts.bonmin.milp_log_level = 0;
cent_opts.bonmin.nlp_log_level = 0;
cent_opts.bonmin.oa_log_level = 0;
cent_opts.print_time = 0;
cent_opts.discrete = discrete_cent;

nlp =   struct('x',x_glob,'f',obj_cent(x_glob),'g',cons_cent(x_glob));
cas =   nlpsol('solver','bonmin',nlp,cent_opts);
sol =   cas('lbx', vertcat(lb{:}),...
            'ubx', vertcat(ub{:}),...
            'lbg', lbg,...
            'ubg', ubg,...
            'x0', vertcat(xz0{:}));
Cent_time = toc;

cas.stats.return_status
Cent_sol = full(sol.x);

%% PaDOA solution
padoa_opts.MINLP_solver = 'cplex';
padoa_opts.MILP_solver = 'cplex';
padoa_opts.MINLPsolver_opts.print_time = 0;
padoa_opts.MINLPsolver_opts.error_on_fail = 0;
padoa_opts.milp_opts.cplex.CPX_PARAM_HEURFREQ = -1;
padoa_opts.milp_opts.error_on_fail = 0;

params.update_hyper_planes_with_U = 1; %options are: 1, 0 
params.remove_inac_hyperplanes = 1; %options are: 1, 0
params.limited_hyperplane_gen = 0; %options are: 1, 0
params.upper_bound_MILP_epigraph = 0;  %options are: 1, 0
params.lower_bound_MILP_epigraph = 1;  %options are: 1, 0

[ PaDOA_xopt, logg ] = PaDOA(obj_funs,cons_funs,lb,ub,A,B,c,x0,z0,params,padoa_opts);

PaDOA_time = sum(logg.time.subtime) + sum(logg.time.milptime) + sum(logg.time.hyptime);

%% display solution
disp('Bonmin time')
disp(Cent_time)
disp('PaDOA time')
disp(PaDOA_time)
disp('PaDOA obj. / Bonmin obj.')
disp(full(obj_cent(PaDOA_xopt)/sol.f))
disp('|PaDOA x - Bonmin x|')
disp(sum(abs(PaDOA_xopt-full(sol.x))))
