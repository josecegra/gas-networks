function [ J ] = jacob_x_b( S )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
v2struct(S)

F0 = F_boundary( x,Nv,FN,sN,Nodes,T_in,p_in,q_out,pipes,w_pipes,e_pipes,alg_pipes,tau,j_t,t_ev,q_in,gasfile);

%Initialize Jacobian
J = zeros(size(F0,1),size(x,1));
per = 1e-5;
for j = 1:size(x,1)
    xt = x;
    dx = per*abs(xt(j));
    xt(j) = xt(j) + dx;
    Ft = F_boundary( xt,Nv,FN,sN,Nodes,T_in,p_in,q_out,pipes,w_pipes,e_pipes,alg_pipes,tau,j_t,t_ev,q_in,gasfile);
    J(:,j) = (Ft - F0)/dx;   
end

J = sparse(J);

end

