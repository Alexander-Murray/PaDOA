close all;
clear all;
clc;

import casadi.*

rand('seed', 0);
randn('seed', 0);

n = 4; %number of features
m = 10; %number of training examples
N = 10; %nuber of subsystems

if mod(m,N)~=0
   error('m must be divisible by N!') 
end

% algorithm parameters
eps = 10^-1; % Termination threshold
max_iter = 500; % Iteration limit
max_time = 3000; % time limit;

w = sprandn(n, 1, 100/n);       % N(0,1), 10% sparse
v = randn(1);                  % random intercept

X0 = sprandn(m*N, n, 10/n);           % data / observations
btrue = sign(X0*w + v);

% noise is function of problem size use 0.1 for large problem
b0 = sign(X0*w + v + sqrt(0.1)*randn(m*N, 1)); % labels with noise

% packs all observations in to an m*N x n matrix
A0 = spdiags(b0, 0, m*N, m*N) * X0;

x_true = [v; w];

% objective and constraints
ratio = sum(b0 == 1)/(m*N);
mu = 0.1*1/(m*N) * norm((1-ratio)*sum(A0(b0==1,:),1) + ratio*sum(A0(b0==-1,:),1), 'inf');
C = [-b0 -A0]';
x = SX.sym('x',[n+1,1]);
% A_cons = diag(ones((n+1)*N,1))-diag(ones((n+1)*N-1,1),1); A_cons(end,1)=-1;
A_cons = eye((n+1)*N) - diag(ones((n+1)*(N-1),1),n+1); A_cons((n+1)*(N-1)+1:(n+1)*N,1:n+1)=-eye(n+1);
obj_cent = 0;
x_glob = SX.sym('x',[(n+1)*N,1]);
% SBB_obj = 0;
for part = 1:N
    c_loc = C(:,1+(part-1)*m:part*m)';
    obj{part} = sum(exp(c_loc*x))+m*mu*norm(x(2:end),1);
    obj_fun{part} = Function('f',{x},{obj{part}});
    lbx{part} = -3*ones(n,1); ubx{part}=3*ones(n,1);
    lbz{part} = -3*ones(1,1); ubz{part}=3*ones(1,1);
    lb{part} = [lbx{part};lbz{part}];
    ub{part} = [ubx{part};ubz{part}];
    loc_cons = [x-ub{part};lb{part}-x];
    con_fun{part} = Function('g',{x},{loc_cons});
    A_con{part} = A_cons(:,(part-1)*(n+1)+1:part*(n+1)-1); %consensus in real vars
    B_con{part} = A_cons(:,part*(n+1)); %consensus in discrete vars
    AB{part} = [A_con{part},B_con{part}];
    discrete{part} = zeros(1,n+1); discrete{part}(n+1)=1;
    x0{part} = zeros(n,1);
    z0{part} = 0;
    u0{part} = [x0{part};z0{part}];
    
    x_loc = x_glob((part-1)*(n+1)+1:part*(n+1));
    obj_cent = obj_cent + obj_fun{part}(x_loc);
%     SBB_obj = SBB_obj + obj_fun{part}([x_glob((part-1)*n+1:part*n);x_glob(n*N+part)]);
end
b_con = zeros((n+1)*N,1);

lbz_glob = vertcat(lbz{:});
ubz_glob = vertcat(ubz{:});

int_parts = [];
p_count = 1;
v_count = 1;
for i = 1:length(lbz_glob)
   ints{i} = lbz_glob(i):ubz_glob(i); 
   int_parts = [int_parts;[i,p_count,v_count]]; %ints{i} belongs to partition p_count and is the v_count^th integer of that partition
   if i==length(vertcat(lbz{1:p_count}))
      p_count = p_count + 1;
      v_count = 1;
   else
       v_count = v_count + 1;
   end
end

%% Solve problem
obj_fun_cent = Function('f',{x_glob},{obj_cent});
for i = 1:N
   AB{i} = [A_con{i} B_con{i}]; 
end
cons_cent = [horzcat(AB{:})*x_glob-b_con; -horzcat(AB{:})*x_glob+b_con];
cons_fun_cent = Function('g',{x_glob},{cons_cent});
opts.discrete = horzcat(discrete{:});
opts.print_time = 0;

tic
nlp =   struct('x',x_glob,'f',obj_cent,'g',cons_cent);
cas =   nlpsol('solver','bonmin',nlp,opts);
sol =   cas('lbx', vertcat(lb{:}),...
            'ubx', vertcat(ub{:}),...
            'lbg', -inf(length(cons_cent),1),...
            'ubg', zeros(length(cons_cent),1),...
            'x0', vertcat(u0{:}));
Cent_time = toc;

%%
% PaDOA options
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

[PADOA_sol, logg, PADOA_time, PADOA_obj] = PaDOA(obj_fun,con_fun,lb,ub,A_con,B_con,b_con,x0,z0,params,padoa_opts);

%% display results
Padoa_iter = size(logg.x_plus,2);
hold on
plot(1:Padoa_iter,full(logg.U),'--','Color',[0 0.4470 0.7410],'LineWidth',3);
plot(1:Padoa_iter,full(logg.L),'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',3);
plot(1:Padoa_iter,full(logg.U),'.','MarkerSize',30,'MarkerEdgeColor',[0 0.4470 0.7410],'MarkerFaceColor','blue');
plot(1:Padoa_iter,full(logg.L),'.','MarkerSize',30,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor','red');
hold off
ylabel('Objective Value')
ylim([full(logg.L(1)-0.1*abs(logg.L(1))),logg.U(1)+0.1*abs(logg.U(1))]);
xlim([1,max(Padoa_iter,2)]);
xlabel('Iterations')
set(gcf,'position',[0 0 1000 420])
set(gca,'fontsize', 14)
box on

disp('Bonmin time')
disp(Cent_time)
disp('PaDOA time')
disp(PADOA_time)
disp('PaDOA iterations')
disp(Padoa_iter)
disp('PaDOA obj. / Bonmin obj.')
disp(full(PADOA_obj/sol.f))