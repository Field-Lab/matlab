
function [C,U,V,stim,cvx_optval]=cvx_polar_admm_sqrt_prox(data,data_event,radii_vec,theta_vec,time_sample_len,dict,lam,vec_len,wave_sample_len,stim_avg,mu,y_lag,gam,A_stim,C_est,prox_lam)
cvx_begin quiet

variables C(vec_len) U(vec_len) V(vec_len) stim(wave_sample_len)
subject to
% 
% y=0*data';
% 
% stim_full=[stim;zeros(time_sample_len,1)];
% 
% for itime=1:time_sample_len
% %     for icell=2:no_cells+1
% %         y = y + circshift(waveforms_full(:,icell),itime)*cell_events(icell-1,itime);
% %     end
% %     
% p=circshift(stim_full,itime-offset);
%     y =  y + p(1:time_sample_len)*data_event(itime); %there is an error here! - look into it! 
% end

y =  dict*vec([C, U, V]') + A_stim*stim;

 norms([U V], 2, 2) <= radii_vec .* C;

    % Linear constraint
    U >= radii_vec .* cos(theta_vec) .* C;

minimize (gam*norm((y-data'),2)  + norm(lam.*C,1)+ (1/(2*prox_lam))*(C-C_est)'*(C-C_est))%+(1/(2*mu))*norm((stim-stim_avg),2))
%minimize (gam*sum(huber(y-data'))  + norm(lam.*C,1)+y_lag'*(stim-stim_avg))%+(1/(2*mu))*norm((stim-stim_avg),2)) 
% NISHAL TODO
 
cvx_end

am_term = (gam*norm((y-data'),2) )
lam_term = norm(lam.*C,1)
y_term = y_lag'*(stim-stim_avg)
mu_term = (1/(2*mu))*norm((stim-stim_avg),2)
