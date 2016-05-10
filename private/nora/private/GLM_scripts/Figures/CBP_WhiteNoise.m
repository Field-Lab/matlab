%% LOAD DATA
%{
load('/Volumes/Analysis/nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_noCP_CBP_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_1772.mat')
fittedGLM_CBP=fittedGLM;
load('/Volumes/Analysis/nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_1772.mat');
%}


load('/Volumes/Analysis/nora/NSEM/GLM_Output/rk1_MU_PS_noCP_CBP_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/OFFPar_1471.mat')
fittedGLM_CBP=fittedGLM;
load('/Volumes/Analysis/nora/NSEM/GLM_Output/rk1_MU_PS_noCP_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/OFFPar_1471.mat');

figure; 

dt = fittedGLM.t_bin;
bins     = 120 * 2 * fittedGLM.bins_per_frame;
time     = dt*[1:bins];

subplot(2,1,1);
rec_rast = fittedGLM.xvalperformance.rasters.recorded(:,1:bins);
rec_rast_CBP = fittedGLM_CBP.xvalperformance.rasters.recorded(:,1:bins);

trials   = size(rec_rast,1);
runawaytrial=ones(trials,1);
ylim([1 trials*2])

hold on
for i_trial = 1:trials
    sim1 = time(find(rec_rast_CBP(i_trial,:)));
    rec1 = time(find(rec_rast(i_trial,:)));
    plot(rec1, i_trial, 'r.')
    plot(sim1, i_trial + trials, 'b.')
end

runtrials=find(runawaytrial);
N=length(runtrials);

convolve=100;
PSTH_rec=zeros(length(runtrials),bins);
PSTH_rec_CBP=zeros(length(runtrials),bins);

for i=1:length(runtrials)
    i_trial=runtrials(i);
    PSTH_rec(i,:)=conv(fittedGLM.xvalperformance.rasters.recorded(i_trial,1:bins),ones(convolve,1),'same');
    PSTH_rec_CBP(i,:)=conv(fittedGLM_CBP.xvalperformance.rasters.recorded(i_trial,1:bins),ones(convolve,1),'same');
end

subplot(2,1,2)
hold on
plot(time,mean(PSTH_rec),'Color',[1 0 0]);
plot(time,mean(PSTH_rec_CBP),'Color',[0 0 1]);
xlim([time(1) time(end)])
hold off 