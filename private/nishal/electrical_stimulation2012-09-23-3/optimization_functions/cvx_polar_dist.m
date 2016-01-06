
function [C,U,V,stim]=cvx_polar_dist(data,data_event,radii_vec,theta_vec,time_sample_len,offset,dict,lam,vec_len,wave_sample_len,stim_avg,mu)
cvx_begin quiet

variables C(vec_len) U(vec_len) V(vec_len) stim(wave_sample_len)
subject to

y=0*data';

stim_full=[stim;zeros(time_sample_len,1)];

for itime=1:time_sample_len
%     for icell=2:no_cells+1
%         y = y + circshift(waveforms_full(:,icell),itime)*cell_events(icell-1,itime);
%     end
%     
p=circshift(stim_full,itime-offset);
    y =  y + p(1:time_sample_len)*data_event(itime); %there is an error here! - look into it! 
end

y = y + (dict)*[C;U;V];

 norms([U V], 2, 2) <= radii_vec .* C;

    % Linear constraint
    U >= radii_vec .* cos(theta_vec) .* C;

minimize (norm(y-data',2)  +lam* norm(C,1)+(1/(2*mu))*sum((stim-stim_avg).^2))

cvx_end

C