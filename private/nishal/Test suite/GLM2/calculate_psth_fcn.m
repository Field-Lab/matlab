
function [PSTH_rec,time]=calculate_psth_fcn(convolve,fittedGLM,rec_rast)

dt = fittedGLM.t_bin;
bins     = fittedGLM.GLMPars.approx_refresh_hz *48 * fittedGLM.bins_per_frame;
time     = dt*[1:bins];



%convolve=150;
PSTH_rec=0*rec_rast;
N=size(rec_rast,1);
for i_trial=1:N
PSTH_rec(i_trial,:)=conv(rec_rast(i_trial,:),ones(convolve,1),'same');
end
PSTH_rec=mean(PSTH_rec);
time=time(1:length(PSTH_rec));

% figure;
% plot(time,PSTH_rec);
end