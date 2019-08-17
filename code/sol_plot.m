qf   = x(1:sum(Nv));
pf   = x(sum(Nv)+1:2*sum(Nv));
for j = 1:w_pipes
    q_wp(j) = x(sN(size(sN,1))+j);
    p_wp(2*j-1) = x(sN(size(sN,1))+w_pipes + 2*j -1);
    p_wp(2*j) = x(sN(size(sN,1))+w_pipes + 2*j); 
end


if ismember('alg',models)
    Re_w = abs(q_wp)*D/(A*eta);
    %lambda_w = 1/4*(log10(k/(3.7065*D) - 5.0425./Re_w.*log10((k/D)^1.1098/2.8257 + 5.8506./Re_w.^0.8981))).^(-2);
    lambda_w = fric*ones(size(q_wp,2),1);
end

p_c     = 46.4e5;                         % Pseudo-critical pressure, see http://link.springer.com/chapter/10.1007/978-3-642-21013-6_15/fulltext.html#Tab1
T_c     = 190.71;
alpha    = 0.257/p_c - 0.533*T_c/(p_c*T);

q_p = [];
p_p = [];
N_vv = [];
for i=1:pipes+w_pipes
    if ismember(i,e_pipes)
        j = find(e_pipes==i);
        N_v = L_e(j)/h+1;
        
        q_p = [q_p;qf(sN(j)+1:sN(j+1))];
        p_p = [p_p;pf(sN(j)+1:sN(j+1))];
        N_vv = [N_vv,N_v];
    else
        j = find(alg_pipes==i);
        N_v = L_w(j)/h+1;
        space = h*(0:N_v-1)';
        q_p = [q_p;q_wp(j)*ones(N_v,1)];
        p_p = [p_p;sqrt(p_wp(2*j-1)^2-((R*T*(1+alpha*(p_wp(2*j)+p_wp(2*j-1))/2))/(D*A^2))*lambda_w(j)*space*q_wp(j)^2)];
        N_vv = [N_vv,N_v];   
    end
end


[hAx,hLine1,hLine2] = plotyy(h*(0:(sum(N_vv)-1)),q_p,h*(0:(sum(N_vv)-1)),p_p);
ylabel(hAx(1),'q') % left y-axis 
ylabel(hAx(2),'p') % right y-axis

