
% Load some data
load('/Volumes/Analysis/nora/NSEM/GLM_Output/rk1_MU_PS_CP_p8IDp8/standardparams/NSEM_mapPRJ/2013-08-19-6/ONPar_2824.mat');
fittedGLM_rk1=fittedGLM;
load('/Volumes/Analysis/nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/NSEM_mapPRJ/2013-08-19-6/ONPar_2824.mat');
fittedGLM_STA=fittedGLM;

% Performance and plotting

dt = fittedGLM_rk1.t_bin;
bins     = 120 * 8 * fittedGLM_rk1.bins_per_frame;
time     = dt*[1:bins];
rec_rast = fittedGLM_rk1.xvalperformance.rasters.recorded(:,1:bins);
sim_rast = fittedGLM_rk1.xvalperformance.rasters.glm_sim(:,1:bins);
trials   = size(rec_rast,1);
    
hFig=figure;
set(hFig, 'Position', [100 100 1000 250])

hold on
runawaytrial=ones(trials,1);
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    sim1 = time(find(sim_rast(i_trial,:)));
    if length(sim1) > 4*length(rec1)
        runawaytrial(i_trial)=0;
    end
end

runtrials=find(runawaytrial);
N=length(runtrials);
convolve=150;
PSTH_rec=zeros(length(runtrials),bins);
PSTH_sim=zeros(length(runtrials),bins);

for i=1:length(runtrials)
    i_trial=runtrials(i);
    PSTH_rec(i,:)=conv(rec_rast(i_trial,:),ones(convolve,1),'same');
    PSTH_sim(i,:)=conv(sim_rast(i_trial,:),ones(convolve,1),'same');
end

hold on
plot(time,mean(PSTH_sim),'r');
xlim([time(1) time(end)])

dt = fittedGLM_STA.t_bin;
bins     = 120 * 8 * fittedGLM_rk1.bins_per_frame;
time     = dt*[1:bins];
rec_rast = fittedGLM_STA.xvalperformance.rasters.recorded(:,1:bins);
sim_rast = fittedGLM_STA.xvalperformance.rasters.glm_sim(:,1:bins);
trials   = size(rec_rast,1);

hold on
runawaytrial=ones(trials,1);
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    sim1 = time(find(sim_rast(i_trial,:)));
    if length(sim1) > 4*length(rec1)
        runawaytrial(i_trial)=0;
    end
end

runtrials=find(runawaytrial);
N=length(runtrials);
convolve=150;
PSTH_rec=zeros(length(runtrials),bins);
PSTH_sim=zeros(length(runtrials),bins);

for i=1:length(runtrials)
    i_trial=runtrials(i);
    PSTH_rec(i,:)=conv(rec_rast(i_trial,:),ones(convolve,1),'same');
    PSTH_sim(i,:)=conv(sim_rast(i_trial,:),ones(convolve,1),'same');
end

hold on
plot(time,mean(PSTH_sim),'b');
xlim([time(1) time(end)])

plot(time,mean(PSTH_rec),'k');

legend('Fit Spatial','STA Spatial','Recorded')