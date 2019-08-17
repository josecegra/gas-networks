x_0_t = x;
%timestep
tau_v = 0.5;
tau = tau_v;

sol_plot

q0   = q_p;
p0   = p_p;

%Number of timesteps
N_t = 50;

%Number of measurements along the network
N_meas = 100;
p_meas = zeros(N_meas,N_t);
q_meas = zeros(N_meas,N_t);

%maximum number of iterations for Newton's method
Ntmax = 300;

%newton method
j_t      = 0;
x_t      = x_0_t;                                 % set starting value
x_t      = double(x_t);
tau_r_s= 1e-4;      % accuracy steady state solution
rel_dif =1;         % Relative difference between x_0 and x
tau_r  = 1e-4;                                % relative error tolerance
t_ev = 1;

figure
time_sol_plot
drawnow
for j_t = 1:N_t
    x_0_t         = x_t;
    vars = {tau,fric,D,eta,k,T,R,p_c,T_c,hprime,g,k_w,T_w,c_v,c,p_En,T_En,q_Ex,x_0_t,x_t};
    F_eval = F(vars{:});
    %Coupling conditions
    F_b = F_boundary( x_t,Nv,FN,sN,Nodes,T_in,p_in,q_out,pipes,w_pipes,e_pipes,alg_pipes,tau,j_t,t_ev,q_in,gasfile);   
    F_eval = [F_eval; F_b];
    r_0    = norm(F_eval);
    i=0;
    while norm(F_eval,inf) >= tau_r * r_0 && i<Ntmax
        J_x_eval    = J_x(vars);
        %Coupling conditions
        S = v2struct;
        J_x_b = jacob_x_b(S);
        J_x_eval = [J_x_eval; J_x_b];
        
        Dx          = J_x_eval\-F_eval;
        x_t           = Dx + x_t;
        vars{1,end} = x_t;
        F_eval      = F(vars{:});
        %Coupling conditions
        F_b = F_boundary( x_t,Nv,FN,sN,Nodes,T_in,p_in,q_out,pipes,w_pipes,e_pipes,alg_pipes,tau,j_t,t_ev,q_in,gasfile);
        F_eval = [F_eval; F_b];
        i           = i+1;
    end
    if any(x_t(sum(Nv)+1:2*sum(Nv))<0)
        time_sol_plot
        error('p became smaller than 0 at iteration j=%d,i=%d at xindex=%d',j_t,i,find(x_t(sum(Nv)+1:2*sum(Nv))<0,1))
    end
    rel_dif=norm((x_t-x_0_t)./x_0_t,inf);
    j_t=j_t+1
    time_sol_plot
    title(j_t)
    drawnow
    
    qf = q_p;
    pf = p_p;
    
    for i = 1:N_meas-1
        p_meas(i,j_t-1) = pf(round((i+1)*size(p_p,1)/N_meas)) - p0(round((i+1)*size(p_p,1)/N_meas));
        q_meas(i,j_t-1) = qf(round((i+1)*size(q_p,1)/N_meas)) - q0(round((i+1)*size(q_p,1)/N_meas));
    end
end

ampl_p = zeros(pipes+w_pipes,1);
ampl_q = zeros(pipes+w_pipes,1);

for i=1:N_meas-1
    ampl_p(i) = max(p_meas(i,:));
    ampl_q(i) = max(q_meas(i,:));
end

plot_perturbation






