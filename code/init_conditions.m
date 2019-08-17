
p = zeros(sum(Nv),1);
q = zeros(sum(Nv),1);
p_alg = [];
q_alg = [];
p_entry = p_in(1);
q_s = q_start(Nodes,pipes+w_pipes,q_out,q_in); % For q we want to fulfill the massflow in the nodes

for i=1:pipes+w_pipes
    if ismember(i,e_pipes)
        j = find(e_pipes==i);
        N_v = L_e(j)/h+1;
        space = h*(0:N_v-1)';
        icp=sN(j)+1:sN(j+1);
        q(icp)        = q_s(j);
        p(icp) = sqrt(p_entry^2-(R*T)/(D*A^2)*fric_v*space*q_s(1)^2);
        p_entry = p(sN(j+1));
    else
        j = find(alg_pipes==i);
        q_alg(j,1) = q_s(pipes+j);
        p_alg(2*j-1,1) = p_entry;
        p_alg(2*j,1) = sqrt(p_entry^2-(R*T)/(D*A^2)*fric_v*L_w(j)*q_s(1)^2);
        p_entry = p_alg(2*j,1);
    end
end