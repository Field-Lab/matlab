function [K_xkp1,B_xkp1] = update_Kj_ASM_population_admm (K_zk,B_zk,u_Kk,u_Bk,rho,X,Y,alpha,dt,scale_cell)

% initialize using thresholded values
K = K_zk; B=B_zk; 

% use fminunc ..
T = size(Y,2);
sta = zeros(size(X,1),1);
Nc = size(Y,1);
b_coeff = gpuArray(zeros(Nc,1));
for icell=1:Nc
 sta = sta + X*(alpha(:,icell)'.*Y(icell,:))'*(T/scale_cell(icell))/T;
 b_coeff(icell) = alpha(:,icell)'*Y(icell,:)'*(T/scale_cell(icell))/T;
end


Kidx =1:length(K_zk);
Bidx = Kidx(end)+1:Kidx(end) + length(B_zk);
x0 = (zeros(length(K_zk) + length(B_zk),1));
x0(Kidx) = gather(K_zk);
x0(Bidx) = gather(B_zk);

optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','off',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on');

[ xfinal,fval ] = fminunc (@(x) optmizeK_admm(x,X,sta,b_coeff,rho,K_zk,u_Kk,Kidx,Bidx,dt,scale_cell),x0,optim_struct);

K_xkp1 = xfinal(Kidx);
B_xkp1 = xfinal(Bidx)';

end