clear
q_in = 0; %keep as 0; 0: pressure as input BC, 1: flow rate as input BC (TO DO)

%friction factor. Chen's friction can be used too below
fric_v = 0.01;
fric = sym('fric');

%timestep
tau_v = 5;
tau      = sym('tau'); 
Ntmax=100; %Maximum number of timesteps 

%Network
%gasfile = '../Networks/one_pipe.mat';
%gasfile = '../Networks/three_pipes.mat';
gasfile = '../Networks/Y-shaped.mat';
%gasfile = '../Networks/circular.mat';

% Load network setup
Nodes=[];
lvars={'L','Nodes','p_in','q_out','T_in','discr'};
load(gasfile,lvars{:});
display(discr)
clear lvars

%models of each pipe
if strcmp('../Networks/one_pipe.mat',gasfile)
    models = {'iso'}; %'iso' or 'semi'
elseif strcmp('../Networks/three_pipes.mat',gasfile)
    models = {'semi','alg','semi'}; %'iso','semi' or 'alg'
elseif strcmp('../Networks/Y-shaped.mat',gasfile)
    models = {'semi','alg','semi'}; %'iso','semi' or 'alg'
elseif strcmp('../Networks/circular.mat',gasfile)
    models = {'semi','alg','semi','alg','semi','semi'}; %'iso','semi' or 'alg'
end


e_pipes = []; %pipes described by transient models
alg_pipes = []; %pipes described by the algebraic model
for i = 1:size(models,2)
    if strcmp(models(i),'iso')||strcmp(models(i),'semi')
        e_pipes = [e_pipes,i];
        
    elseif strcmp(models(i),'alg')
        alg_pipes = [alg_pipes,i];
    end
end
L_e = L(e_pipes); %Length of the transient pipes
L_w = L(alg_pipes); %length of the algebraic pipes


h = 100; %lattice spacing
validateattributes(L_e/h,{'double'},{'integer','positive'}) %check if L and h produce integers
pipes = length(L_e);
Ni=L_e/h; % Number of intervals per pipe
Nv=L_e/h+1; % Number of variables per pipe
w_pipes = length(L_w);
numVars=2;

%indexing
sN = zeros(numVars*pipes+1,1);    % Vector with the indices were the start of q,p,rho,F,T per pipe can be found
for j=1:pipes
    for i=0:numVars-1 % Indices over the number of calculated output variables (q,p,rho,F,T)
        sN(pipes*i+j) = i*sum(Nv)+sum(Nv(1:j-1));
    end
end
sN(end) = numVars*sum(Nv); % Last index of solutions
FN=zeros(numVars*pipes+2,1);
for i=0:1
    for j=1:pipes
        FN(pipes*i+j+1)= FN(pipes*i+j)+Ni(j);
    end
end
FN(end)=numVars*sum(Nv); %last equation


%declare symbolic
% Solution and independent/uncertain variables
x      = sym('x', [sN(end)+3*w_pipes 1],'real');
x_0    = sym('x_0', [sN(end)+3*w_pipes 1],'real');
c      = sym('c');                      % Speed of sound
p_En   = sym('p_En', [pipes 1]);                   % Entrance pressure
q_Ex   = sym('q_Ex', [pipes 1]);            % Exit q
T_En   = sym('T_En', [pipes 1]);                   % Entrance temperature
D      = sym('D');
k      = sym('k');                      % Roughness
eta    = sym('eta');                    % Dynamic viscosity
T      = sym('T');                                 % Constant temperature
R      = sym('R');                      % Specific gas constant
p_c    = sym('p_c');                    % Pseudo-critical pressure
T_c    = sym('T_c');                    % Pseudo-critical temperature
hprime = sym('hprime');                 % Slope of the pipe
g      = sym('g');                      % Acceleration due to gravity
k_w    = sym('k_w');                    % Heat conductivity coefficient
T_w    = sym('T_w');                    % Pipe wall temperature
c_v    = sym('c_v');                    % Volumetric heat capacity

% Vectors with multiple values for pipe dependent variables
Dv=ones(pipes,1)*D;
kv=ones(pipes,1)*k;
hprimev=ones(pipes,1)*hprime;
T_wv=ones(pipes,1)*T_w;
hv=ones(pipes,1)*h;
% Dependent variables
A      = pi*(D/2).^2;                    % Cross-sectional area [m^2]
Av      = pi*(Dv/2).^2;                    % Cross-sectional area [m^2]


iq_1   = (1:sum(Nv))';      % Indices in x for q_{j,i-1}
iq_1(sN(2:pipes+1))=[];
iq     = iq_1+1;  % Indices in x for q_{j,i}
ip     = iq+sum(Nv);    % Indices in x for p_{j,i}
ip_1   = iq_1+sum(Nv);    % Indices in x for p_{j,i-1}

p_0   = sym('p_0', [sum(Nv)-pipes 1]);        % Used to set the size of p_0 and all other initial values
p_0_1 = p_0;  p       = p_0;  p_1   = p_0;
q_0   = p_0;  q_0_1   = p_0;  q     = p_0;  q_1   = p_0;

Re    = p_0;  Re_1=p_0; lambda=p_0; lambda_1=p_0;

for j=1:pipes
    icp            = sN(j)+2-j:sN(j+1)-j;       % Indices of the current pipe
    
    q_0_1(icp)     = x_0(iq_1(icp));            % q^0_{j,i-1}
    q_1(icp)       = x(iq_1(icp));              % q^1_{j,i-1}
    p_0(icp)       = x_0(ip(icp));              % p^0_{j,i}
    p(icp)         = x(ip(icp));                % p^1_{j,i}
    q_0(icp)       = x_0(iq(icp));                     % q^0_{1,i}
    q(icp)         = x(iq(icp));                       % q^1_{1,i}
    p_0_1(icp)     = x_0(ip_1(icp));                   % p^0_{j,i-1}
    p_1(icp)       = x(ip_1(icp));                     % p^1_{j,i-1}
    
    Re(icp)       = abs(q(icp))*Dv(j)/(Av(j)*eta); % Reynolds number with q^1_{j,i}
    Re_1(icp)     = abs(q_1(icp))*Dv(j)/(Av(j)*eta);              % Reynolds number with q^1_{j,i-1}
    
    % friction factor
    %lambda(icp)   = 1/4*(log10(kv(j)/(3.7065*Dv(j)) - 5.0425./Re(icp).*log10((kv(j)/Dv(j))^1.1098/2.8257 + 5.8506./Re(icp).^0.8981))).^(-2);
    lambda(icp)   = fric;

    % lambda with Re_1
    %lambda_1(icp) = 1/4*(log10(kv(j)/(3.7065*Dv(j)) - 5.0425./Re_1(icp).*log10((kv(j)/Dv(j))^1.1098/2.8257 + 5.8506./Re_1(icp).^0.8981))).^(-2);
    lambda_1(icp) = fric;
      
end

% Define the full variable vectors
qf   = x(1:sum(Nv));
pf   = x(sum(Nv)+1:2*sum(Nv));
q_w   = sym('q_w', [w_pipes 1]); 
p_w   = sym('p_w', [2*w_pipes 1]); 
for j = 1:w_pipes
    q_w(j) = x(sN(end)+j);
    p_w(2*j-1) = x(sN(end)+w_pipes + 2*j -1);
    p_w(2*j) = x(sN(end)+w_pipes + 2*j); 
    
end

% Assemble F
F = sym('F', [FN(end-1) 1]);

alpha    = 0.257/p_c - 0.533*T_c/(p_c*T);   % Used in z(p) = 1 + alpha*p
Cv      = R*T./Av;

%set friction factor algebraic model
if ismember('alg',models)
    Re_w = abs(q_w)*D/(A*eta);
    %lambda_w = 1/4*(log10(k/(3.7065*D) - 5.0425./Re_w.*log10((k/D)^1.1098/2.8257 + 5.8506./Re_w.^0.8981))).^(-2);
    lambda_w = fric*ones(size(q_w,1),1);
end

%build nonlinear equations for transient models
for i=1:pipes+w_pipes
    if strcmp(models(i),'iso')||strcmp(models(i),'semi')
        j = find(e_pipes==i);
        icp = FN(j)+1:FN(j+1);  % Indices of the current pipe
        
        if strcmp(models(i),'iso')       
            % Continuity equation
            F(FN(j)+1:FN(j+1))   = (p_1(icp)./(1+alpha*p_1(icp)) ...
                + p(icp)./(1+alpha*p(icp)))/(2*tau) ...
                - (p_0_1(icp)./(1+alpha*p_0_1(icp)) ...
                + p_0(icp)./(1+alpha*p_0(icp)))/(2*tau) ...
                + Cv(j)*(q(icp) - q_1(icp))/hv(j);
            
            % Momentum equation
            F(FN(j+pipes)+1:FN(j+pipes+1)) = (q_1(icp) + q(icp))/(2*tau) ...
                - (q_0_1(icp) + q_0(icp))/(2*tau) ...
                + Cv(j)/hv(j)*((1+alpha*p(icp)).*q(icp).^2./p(icp) ...
                - (1+alpha*p_1(icp)).*q_1(icp).^2./p_1(icp)) ...
                + Av(j)*(p(icp) - p_1(icp))/h ...
                + g*hprimev(j)*(p_1(icp)./(1+alpha*p_1(icp)) ...
                + p(icp)./(1+alpha*p(icp)))/(Cv(j)*2) ...
                + Cv(j)/(4*Dv(j))*((1+alpha*p_1(icp)).*q_1(icp).*abs(q_1(icp)).*lambda_1(icp)./p_1(icp) ...
                + (1+alpha*p(icp)).*q(icp).*abs(q(icp)).*lambda(icp)./p(icp));
            
        elseif strcmp(models(i),'semi')
            % Continuity equation
            F(FN(j)+1:FN(j+1))   = (p_1(icp)./(1+alpha*p_1(icp)) ...
                + p(icp)./(1+alpha*p(icp)))/(2*tau) ...
                - (p_0_1(icp)./(1+alpha*p_0_1(icp)) ...
                + p_0(icp)./(1+alpha*p_0(icp)))/(2*tau) ...
                + Cv(j)*(q(icp) - q_1(icp))/hv(j);
            
            % Momentum equation
            F(FN(j+pipes)+1:FN(j+pipes+1)) = (q_1(icp) + q(icp))/(2*tau) ...
                - (q_0_1(icp) + q_0(icp))/(2*tau) ...
                + Av(j)*(p(icp) - p_1(icp))/h ...
                + g*hprimev(j)*(p_1(icp)./(1+alpha*p_1(icp)) ...
                + p(icp)./(1+alpha*p(icp)))/(Cv(j)*2) ...
                + Cv(j)/(4*Dv(j))*((1+alpha*p_1(icp)).*q_1(icp).*abs(q_1(icp)).*lambda_1(icp)./p_1(icp) ...
                + (1+alpha*p(icp)).*q(icp).*abs(q(icp)).*lambda(icp)./p(icp));
                        
        end 
    end
end

%build nonlinear equations for the algebraic model
for j = 1:w_pipes
    F(FN(end-1)+j) = p_w(2*j)^2-p_w(2*j-1)^2 + ((R*T*(1+alpha*(p_w(2*j)+p_w(2*j-1))/2)*L_w(j)*lambda_w(j))/(D*A^2))*q_w(j)^2; 
end

% Internal equations indices
iind=zeros(pipes*2,2);
for j=1:pipes
    iind((j-1)*2+1:j*2,:)=[FN(j)+1, FN(j+1);
        FN(j+1*pipes)+1, FN(j+1*pipes+1)];
end


%Jacob F handles
% Compute the Jacobians with symbolic toolbox
disp('Computing Jacobians w.r.t. x...');
J_xi = jacobian(F, x);
disp('Finished computing the Jacobians.');
vars = {tau,fric,D,eta,k,T,R,p_c,T_c,hprime,g,k_w,T_w,c_v,c,p_En,T_En,q_Ex,x_0,x};
disp('Creating Matlab function F...');
F4t=F; %copy of symbolic F
F = matlabFunction(F,'vars',vars);
disp('Creating Matlab function J_xi...');
J_xi4t = J_xi;
J_xi = matlabFunction(J_xi,'vars',vars);% Jacobian for the internal equations
disp('Finished creating Matlab functions.');

%set values
% Parameter values
fric = fric_v;
tau = tau_v;
p_En    = reshape(p_in,pipes+w_pipes,1);                            % Boundary values
q_Ex    = reshape(q_out,pipes+w_pipes,1);
T_En    = reshape(T_in,pipes+w_pipes,1);
c       = 382.75;                         % Speed of sound [m/s] =sqrt(R*T_0*z_0)=sqrt(518.3*288.15*0.928)
D       = 0.6;                            % Diameter [m] =800 mm
k       = 1.5e-5;                         % Roughness [m]
eta     = 1e-5;                           % Dynamic viscosity [Pa s]
TC       = 288.15;                         % Constant temperature [K]
R       = 500;                          % Specific gas constant [m^2 s^{-2} K^{-1}]
p_c     = 46.4e5;                         % Pseudo-critical pressure, see http://link.springer.com/chapter/10.1007/978-3-642-21013-6_15/fulltext.html#Tab1
T_c     = 190.71;                         % Pseudo-critical temperature, see same table as above
hprime  = 0.0;                           % Slope of the pipe
g       = 9.80665;                        % Acceleration due to gravity [m/s^2]

Dv=ones(pipes,1)*D;
kv=ones(pipes,1)*k;
hprimev=ones(pipes,1)*hprime;
T_wv=ones(pipes,1)*T_w;
hv=ones(pipes,1)*h;

% Dependent variables
A      = pi*(D/2).^2;                    % Cross-sectional area [m^2]
Av      = pi*(Dv/2).^2;                    % Cross-sectional area [m^2]
Cv      = R*TC./Av;
clear TC

%Temperature
T = 293;
dummy = 0;
init_conditions %set initial conditions
%initial state
x_0 = [q;p;q_alg;p_alg];
%build the jacobian
J_x=@(vars) ...
    J_builder(J_xi,iind,pipes,numVars,vars,FN);

%newton method
j_t      = 0;
x      = x_0;                                 % set starting value
x      = double(x);
tau_r_s= 1e-3;      % accuracy steady state solution
rel_dif =1;         % Relative difference between x_0 and x
tau_r  = 1e-4; % relative error tolerance
t_ev = 0;

figure
sol_plot
drawnow
while j_t<Ntmax && rel_dif > tau_r_s
    x_0         = x;
    [p_p_0,q_p_0] = get_p_q(x_0,h,Nv,pipes,w_pipes,sN,models,D,A,R,T,eta,k,e_pipes,alg_pipes,L_e,L_w);
    vars = {tau,fric,D,eta,k,T,R,p_c,T_c,hprime,g,k_w,T_w,c_v,c,p_En,T_En,q_Ex,x_0,x};
    F_eval = F(vars{:});
    %Coupling conditions
    F_b = F_boundary( x,Nv,FN,sN,Nodes,T_in,p_in,q_out,pipes,w_pipes,e_pipes,alg_pipes,tau,j_t,t_ev,q_in,gasfile);   
    F_eval = [F_eval; F_b];
    r_0    = norm(F_eval);
    i=0;
    while norm(F_eval) >= tau_r * r_0 && i<Ntmax/10
        J_x_eval    = J_x(vars);
        %Coupling conditions
        S = v2struct;
        J_x_b = jacob_x_b(S);
        J_x_eval = [J_x_eval; J_x_b];
        
        Dx          = J_x_eval\-F_eval;
        x           = Dx + x;
        vars{1,end} = x;
        F_eval      = F(vars{:});
        %Coupling conditions
        F_b = F_boundary( x,Nv,FN,sN,Nodes,T_in,p_in,q_out,pipes,w_pipes,e_pipes,alg_pipes,tau,j_t,t_ev,q_in,gasfile);
        F_eval = [F_eval; F_b];
        i           = i+1; 
    end
    if any(x(sum(Nv)+1:2*sum(Nv))<0)
        sol_plot
        error('p became smaller than 0 at iteration j=%d,i=%d at xindex=%d',j,i,find(x(sum(Nv)+1:2*sum(Nv))<0,1))
    end
    j_t=j_t+1
    sol_plot
    title(j_t)
    drawnow
    [p_p,q_p] = get_p_q(x,h,Nv,pipes,w_pipes,sN,models,D,A,R,T,eta,k,e_pipes,alg_pipes,L_e,L_w);
    rel_dif_p=norm((p_p-p_p_0)./p_p_0,inf);
    rel_dif_q=norm((q_p-q_p_0)./q_p_0,inf);
    rel_dif = max(rel_dif_p,rel_dif_q);
end

if j_t==Ntmax
    warning('Stopped after %d iterations',j_t)
end
fprintf('Total number of steps: %d\n',j_t);
fprintf('Total number of Newton iterations: %d\n',i);
%export the results
qf   = x(1:sum(Nv));
pf   = x(sum(Nv)+1:2*sum(Nv));







