function [fval,grad] = minL_dec(x,K,B,Y,rho,Z_k,u_k,dt,T,scale_cells)


%% objective value

SU_inp = K'*x;
fval = 0;
for icell = 1:size(Y,1)
    cell_su_inp = exp(SU_inp + B(:,icell));
fval = fval+(sum(cell_su_inp)*dt/T - Y(icell,:)*log(sum(cell_su_inp))/T)*(T/scale_cells(icell));
end
fval = fval + (rho/2)*norm(x - Z_k + u_k,2)^2;

%% grad

grad = zeros(size(K,1),1);
for icell=1:size(Y,1)
    cell_su_inp = exp(SU_inp + B(:,icell));
    grad = grad + (K*cell_su_inp*dt/T -Y(icell,:)*K*(cell_su_inp/sum(cell_su_inp))/T)*T/scale_cells(icell);
end
grad = grad + (rho)*(x - Z_k + u_k);


grad=gather(grad);
fval= gather(fval);
end