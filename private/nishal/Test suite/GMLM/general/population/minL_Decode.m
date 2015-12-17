function [X_kp1,f_total] = minL_Decode(K,B,Y,rho,Z_k,u_k,X_k,dt,scale_cells)
T = size(Y,2);

optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','off',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off');

X_kp1 = 0*X_k;

K= gather(K);
B=gather(B);
fval=zeros(T,1);
u_k=gather(u_k);Z_k=gather(Z_k);
tic;
parfor itime = 1:T
[X_kp1(:,itime),fval(itime)] = fminunc(@(x)minL_dec(x,K,B,Y(:,itime),rho,Z_k(:,itime),u_k(:,itime),dt,T,scale_cells),X_k(:,itime),optim_struct);
end
toc;
f_total=sum(fval);

end