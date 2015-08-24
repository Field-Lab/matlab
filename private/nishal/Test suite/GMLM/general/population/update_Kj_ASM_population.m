function [fval,Knew] = update_Kj_ASM_population(X,Y,Bs,Ks,T,Nc,SU_inp,alpha,dt)

sta = zeros(size(X,1),1);
for icell=1:Nc
 sta = sta + X*(alpha(icell,:).*Y(icell,:))'/T;
end

const = sum(exp(Bs))*dt/T;

optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','off',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on');

[ Knew,fval ] = fminunc (@(k) optmizeK(k,X,const,sta),Ks,optim_struct);

end