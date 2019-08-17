function [ Fb ] = F_boundary( x,Nv,FN,sN,Nodes,T_in,p_in,q_out,pipes,w_pipes,e_pipes,alg_pipes,tau,j_t,t_ev,q_in,gasfile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

qf   = x(1:sum(Nv));
pf   = x(sum(Nv)+1:2*sum(Nv));

for j = 1:w_pipes
    q_w(j) = x(sN(end)+j);
    p_w(2*j-1) = x(sN(end)+w_pipes + 2*j -1);
    p_w(2*j) = x(sN(end)+w_pipes + 2*j); 
end

c = 382.75;
if t_ev == 1   
    amp = 10^3;
    p_in(1) = 5e6 + amp*exp(-((c*j_t*tau)^2)/(2*500^2)); 
end

Fb = zeros(sN(end)-FN(end-1),1);
eqc = 1;

for m=1:size(Nodes,1)
    ec=~cellfun('isempty',Nodes(m,:))*[1i; 1]; % if real inflow, imaginary outflow, complex node
    switch ec
        case 1
            % Inflow
            pind=Nodes{m,2};
            if length(pind) ~=1; warning('Node row %d contains multiple inflow nodes',m); end
            if ismember(pind,e_pipes)
                i = find(e_pipes==pind);
                Fb(eqc)=pf(sN(i)+1)-p_in(pind);
            else
                i = find(alg_pipes==pind);
                Fb(eqc)=p_w(2*i-1)-p_in(pind);
            end
            eqc=eqc+1;
              
        case 1i
            if q_in == 0
                % Outflow
                pind=Nodes{m,1};
                if length(pind) ~=1; warning('Node row %d contains multiple outflow nodes',m); end
                if ismember(pind,e_pipes)
                    i = find(e_pipes==pind);
                    Fb(eqc) = qf(sN(i+1))-q_out(pind);
                else
                    i = find(alg_pipes==pind);
                    Fb(eqc) = q_w(i)-q_out(pind);
                end
                eqc=eqc+1;
            end
            
        case 1i+1
            % Real node
            pin=Nodes{m,1};
            pout=Nodes{m,2};
            % set q equation
            Fb(eqc)=0;
            for mk=1:length(pin)
                if ismember(pin(mk),e_pipes)
                    i = find(e_pipes==pin(mk));
                    Fb(eqc)=Fb(eqc)+qf(sN(i+1)); %Sum of flows
                else
                    i = find(alg_pipes==pin(mk));
                    Fb(eqc)=Fb(eqc)+q_w(i); %Sum of flows
                end
            end
            
            for mk=1:length(pout)
                if ismember(pout(mk),e_pipes)
                    i = find(e_pipes==pout(mk));
                    Fb(eqc)=Fb(eqc)-qf(sN(i)+1); % Update the sum of flows
                else
                    i = find(alg_pipes==pout(mk));
                    Fb(eqc)=Fb(eqc)-q_w(i); % Update the sum of flows
                end      
            end
            %Reference pressure
            if ismember(pin(1),e_pipes)
                i = find(e_pipes==pin(1));
                ref_p = pf(sN(i+1));
            else
                i = find(alg_pipes==pin(1));
                ref_p = p_w(2*i);
            end
            %Impose all the input pipes the same ending pressure
            for mk=2:length(pin)
                if ismember(pin(mk),e_pipes)
                    i = find(e_pipes==pin(mk));
                    Fb(eqc+mk-1) = ref_p - pf(sN(i+1));
                else
                    i = find(alg_pipes==pin(mk));
                    Fb(eqc+mk-1) = ref_p - p_w(2*i);
                end
            end
            cp=eqc+(length(pin)-1);
            for mk=1:length(pout)
                if ismember(pout(mk),e_pipes)
                    i = find(e_pipes==pout(mk));
                    Fb(cp+mk)=ref_p-pf(sN(i)+1);
                else
                    i = find(alg_pipes==pout(mk));
                    Fb(cp+mk)=ref_p-p_w(2*i-1);
                end
            end
            eqc=cp+length(pout)+1;
    end
end
                    
                
% Internal equations indices
iind=zeros(pipes*2,2);
for j=1:pipes
    iind((j-1)*2+1:j*2,:)=[FN(j)+1, FN(j+1);
        FN(j+1*pipes)+1, FN(j+1*pipes+1)];
end

end

