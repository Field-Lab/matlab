function [filters,bias]= fit_filter_cvx(y_tsp,kx_filter,tsp,kx_sum,N,mov_filtered_b,dt,mov_filtered_spiked,initf,initb)

xx = y_tsp.*kx_filter(tsp)'./kx_sum;
sta_f= mov_filtered_spiked*xx'/N;
sta_f_b = [sta_f;sum(xx')/N];

optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','off',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on');
% Maybe, hessian has made it slow in past!

 [x,fval,exitflag,output,grad,hessian]  = fminunc(@(x)GMLM_decomposed_fcn(x,mov_filtered_b,sta_f_b,N,dt),[initf';initb],optim_struct);
filters = x(1:end-1);
bias=x(end);

end