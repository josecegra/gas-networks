figure

X = 1:size(p_meas,2);
Y = 1:size(p_meas,1);
Z = p_meas;
contourf(Z,20,':')
colorbar
