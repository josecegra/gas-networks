function q_s = q_start(Nodes,pipes,q_out,q_in)
% Calculate the initial values for q for each pipe

q_s=nan(pipes,1);
lc=1;

while any(isnan(q_s)) && lc <= size(Nodes,1)^2
    m=mod(lc-1,size(Nodes,1))+1;
    ec=~cellfun('isempty',Nodes(m,:))*[1i; 1]; % if real inflow, imaginary outflow, complex node
    
    switch ec
        case 1
            % Inflow
            % do nothing
            if q_in ==1
                pind=Nodes{m,2};
                if length(pind) ~=1; warning('Node row %d contains multiple outflow nodes',m); end
                q_s(pind)=q_out(pind);
                
            end
        case 1i
            % Outflow
            if q_in == 0
                pind=Nodes{m,1};
                if length(pind) ~=1; warning('Node row %d contains multiple outflow nodes',m); end
                q_s(pind)=q_out(pind);
            end
            
        case 1i+1
            % Real node
            pin=Nodes{m,1};
            pout=Nodes{m,2};
            
            if any(isnan(q_s(pin)))&& ~any(isnan(q_s(pout)))
                q_s(pin)=sum(q_s(pout))/length(pin);
            end
        otherwise
            warning('Empty node row %d provided',m)
    end

    lc=lc+1;
end
end

