
function [PSTH_rec,time]=calculate_psth_fcn2(convolve,dt,len,rec_rast)


bins     = len;
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