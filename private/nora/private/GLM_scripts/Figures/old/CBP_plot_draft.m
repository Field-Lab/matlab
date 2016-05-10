
%{
load('/Volumes/Analysis/nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_noCP_CBP_p8IDp8/standardparams/NSEM_mapPRJ/2012-08-09-3/ONPar_1772.mat')
fittedGLM_CBP=fittedGLM;
load('/Volumes/Analysis/nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/NSEM_mapPRJ/2012-08-09-3/ONPar_1772.mat')
%}

% compare the filters
% text(-.1, 1-0.1*c,sprintf('Raw xval bitsperspike %1.3e',performance.logprob_glm_bpspike))
dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;

subplot(2,3,1)
hold on
bins    = [1:length(fittedGLM.linearfilters.PostSpike.Filter)];
time_msec = 1000*dt*bins ;
oneline = ones(1,length(time_msec)); 
plot(time_msec, oneline, 'k-'); hold on
plot(time_msec,exp(fittedGLM.linearfilters.PostSpike.Filter),'black--')
plot(time_msec,exp(fittedGLM_CBP.linearfilters.PostSpike.Filter),'black')
  xlim([0, time_msec(end)]);
    ylim([0, max(1.5, max(exp(fittedGLM.linearfilters.PostSpike.Filter)))]);
    ylabel('gain'); xlabel('msec'); title('Post Spike Filter')
    
subplot(2,3,4)
hold on
frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec,fittedGLM.linearfilters.Stimulus.time_rk1,'black--')
plot(time_msec,fittedGLM_CBP.linearfilters.Stimulus.time_rk1,'black')
plot(time_msec, zeroline, 'k-')
xlim([0, time_msec(end)])
xlabel('msec'); title('Time Filter');
set(gca, 'ytick', 0); 
legend('Vision Sorted Spikes','CBP Sorted Spikes')

subplot(2,3,[2,3]);
title('Vision Sorted Rasters')
bins     = 120 * 8 * fittedGLM.bins_per_frame;
rec_rast = fittedGLM.xvalperformance.rasters.recorded(:,1:bins);
sim_rast = fittedGLM.xvalperformance.rasters.glm_sim(:,1:bins); 
trials   = size(rec_rast,1);
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([1 , 2*trials]); hold on

runawaytrial=ones(trials,1);
for i_trial = 1:trials
    sim1 = time(find(sim_rast(i_trial,:)));
    rec1 = time(find(rec_rast(i_trial,:)));
    plot(rec1, i_trial, 'k.')
    if length(sim1) < 4*length(rec1) 
        plot(sim1, i_trial + trials, 'r.')
        runawaytrial(i_trial)=0;
    end
end

subplot(2,3,[5,6]);
title('CBP Sorted Rasters')
bins     = 120 * 8 * fittedGLM.bins_per_frame;
rec_rast = fittedGLM_CBP.xvalperformance.rasters.recorded(:,1:bins);
sim_rast = fittedGLM_CBP.xvalperformance.rasters.glm_sim(:,1:bins); 
trials   = size(rec_rast,1);
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([1 , 2*trials]); hold on

for i_trial = 1:trials
    sim1 = time(find(sim_rast(i_trial,:)));
    rec1 = time(find(rec_rast(i_trial,:)));
    plot(rec1, i_trial, 'k.')
    if length(sim1) < 4*length(rec1) 
        plot(sim1, i_trial + trials, 'r.')
        runawaytrial(i_trial)=0;
    end
end
runtrials=find(runawaytrial);

bins     = 120 * 8 * fittedGLM.bins_per_frame;
convolve=100;
leftover=9600-convolve+1;
PSTH_rec=zeros(length(runtrials),leftover);
PSTH_sim=zeros(length(runtrials),leftover);
PSTH_rec_CBP=zeros(length(runtrials),leftover);
PSTH_sim_CBP=zeros(length(runtrials),leftover);

for i=1:length(runtrials)
    i_trial=runtrials(i)
    PSTH_rec(i,:)=conv(fittedGLM.xvalperformance.rasters.recorded(i_trial,1:bins),ones(convolve,1),'valid');
    PSTH_sim(i,:)=conv(fittedGLM.xvalperformance.rasters.glm_sim(i_trial,1:bins),ones(convolve,1),'valid');
    PSTH_rec_CBP(i,:)=conv(fittedGLM_CBP.xvalperformance.rasters.recorded(i_trial,1:bins),ones(convolve,1),'valid');
    PSTH_sim_CBP(i,:)=conv(fittedGLM_CBP.xvalperformance.rasters.glm_sim(i_trial,1:bins),ones(convolve,1),'valid');
end

figure(2)
title('PSTH')
hold on
time     = dt*[1:bins];
% fill([time(1:leftover) time(leftover:-1:1)],[mean(PSTH_rec)+std(PSTH_rec) flip(mean(PSTH_rec)-std(PSTH_rec))],'k','FaceAlpha',0.5)
% fill([time(1:leftover) time(leftover:-1:1)],[mean(PSTH_sim)+std(PSTH_sim) flip(mean(PSTH_sim)-std(PSTH_sim))],'r','FaceAlpha',0.5)
% fill([time(1:leftover) time(leftover:-1:1)],[mean(PSTH_rec_CBP)+std(PSTH_rec_CBP) flip(mean(PSTH_rec_CBP)-std(PSTH_rec_CBP))],'k','FaceAlpha',0.5)
% fill([time(1:leftover) time(leftover:-1:1)],[mean(PSTH_sim_CBP)+std(PSTH_sim_CBP) flip(mean(PSTH_sim_CBP)-std(PSTH_sim_CBP))],'r','FaceAlpha',0.5)
plot(time(2000:leftover),mean(PSTH_rec(:,2000:leftover))-mean(mean(PSTH_rec(:,2000:leftover))),'black-.','LineWidth',2)
plot(time(2000:leftover),mean(PSTH_sim(:,2000:leftover))-mean(mean(PSTH_sim(:,2000:leftover))),'red-.','LineWidth',2)
plot(time(2000:leftover),mean(PSTH_rec_CBP(:,2000:leftover))-mean(mean(PSTH_rec_CBP(:,2000:leftover))),'black')
plot(time(2000:leftover),mean(PSTH_sim_CBP(:,2000:leftover))-mean(mean(PSTH_sim_CBP(:,2000:leftover))),'red')
legend('Vision spikes','GLM with Vision spikes', 'CBP spikes','GLM with CBP spikes')
mean(mean(PSTH_rec(:,2000:leftover)))
mean(mean(PSTH_sim(:,2000:leftover)))
mean(mean(PSTH_rec_CBP(:,2000:leftover)))
mean(mean(PSTH_sim_CBP(:,2000:leftover)))
hold off

%% LOAD DATA

load('/Volumes/Analysis/nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_1772.mat')
fittedGLM_CBP=fittedGLM;
load('/Users/Nora/Downloads/ONPar_1772.mat');

%% FIGURE 1: FILTERS

dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;

figure(1)

subplot(2,1,1)
hold on
bins    = [1:length(fittedGLM.linearfilters.PostSpike.Filter)];
time_msec = 1000*dt*bins ;
oneline = ones(1,length(time_msec)); 
plot(time_msec, oneline, 'k-'); hold on
plot(time_msec,exp(fittedGLM.linearfilters.PostSpike.Filter),'black--')
plot(time_msec,exp(fittedGLM_CBP.linearfilters.PostSpike.Filter),'black')
  xlim([0, time_msec(end)]);
    ylim([0, max(1.5, max(exp(fittedGLM.linearfilters.PostSpike.Filter)))]);
    ylabel('gain'); xlabel('msec'); title('Post Spike Filter')
    
subplot(2,1,2)
hold on
frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec,fittedGLM.linearfilters.Stimulus.time_rk1,'black--')
plot(time_msec,fittedGLM_CBP.linearfilters.Stimulus.time_rk1,'black')
plot(time_msec, zeroline, 'k-')
xlim([0, time_msec(end)])
xlabel('msec'); title('Time Filter');
set(gca, 'ytick', 0); 
legend('Vision Sorted Spikes','CBP Sorted Spikes')

%% FIGURE 2: RASTERS

dt = fittedGLM.t_bin;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
time     = dt*[1:bins];

subplot(5,1,1);
title('Vision Sorted Rasters')
rec_rast = fittedGLM.xvalperformance.rasters.recorded(:,1:bins);
sim_rast = fittedGLM.xvalperformance.rasters.glm_sim(:,1:bins); 

trials   = size(rec_rast,1);
runawaytrial=ones(trials,1);
ylim([1 trials*2])

hold on
for i_trial = 1:trials
    sim1 = time(find(sim_rast(i_trial,:)));
    rec1 = time(find(rec_rast(i_trial,:)));
    plot(rec1, i_trial, 'k.')
    if length(sim1) < 4*length(rec1) 
        plot(sim1, i_trial + trials, 'r.')
    else
        runawaytrial(i_trial)=0;
    end
end

subplot(5,1,2);
title('CBP Sorted Rasters')
rec_rast = fittedGLM_CBP.xvalperformance.rasters.recorded(:,1:bins);
sim_rast = fittedGLM_CBP.xvalperformance.rasters.glm_sim(:,1:bins); 
ylim([1 trials*2])
hold on
for i_trial = 1:trials
    sim1 = time(find(sim_rast(i_trial,:)));
    rec1 = time(find(rec_rast(i_trial,:)));
    plot(rec1, i_trial, 'k.')
    if length(sim1) < 4*length(rec1) 
        plot(sim1, i_trial + trials, 'r.')
    else
        runawaytrial(i_trial)=0;
    end
end

runtrials=find(runawaytrial);
N=length(runtrials);

convolve=100;
PSTH_rec=zeros(length(runtrials),bins);
PSTH_sim=zeros(length(runtrials),bins);
PSTH_rec_CBP=zeros(length(runtrials),bins);
PSTH_sim_CBP=zeros(length(runtrials),bins);

for i=1:length(runtrials)
    i_trial=runtrials(i);
    PSTH_rec(i,:)=conv(fittedGLM.xvalperformance.rasters.recorded(i_trial,1:bins),ones(convolve,1),'same');
    PSTH_sim(i,:)=conv(fittedGLM.xvalperformance.rasters.glm_sim(i_trial,1:bins),ones(convolve,1),'same');
    PSTH_rec_CBP(i,:)=conv(fittedGLM_CBP.xvalperformance.rasters.recorded(i_trial,1:bins),ones(convolve,1),'same');
    PSTH_sim_CBP(i,:)=conv(fittedGLM_CBP.xvalperformance.rasters.glm_sim(i_trial,1:bins),ones(convolve,1),'same');
end

subplot(5,1,3:5)
title('PSTH')
hold on
a=fill([time flip(time)],[mean(PSTH_rec)+std(PSTH_rec)/sqrt(N) flip(mean(PSTH_rec)-std(PSTH_rec)/sqrt(N))],[1 0 0],'FaceAlpha',0.2,'EdgeColor','none');
fill([time flip(time)],[mean(PSTH_sim)+std(PSTH_sim)/sqrt(N) flip(mean(PSTH_sim)-std(PSTH_sim)/sqrt(N))],[1 0 0],'FaceAlpha',0.5,'EdgeColor','none');
b=fill([time flip(time)],[mean(PSTH_rec_CBP)+std(PSTH_rec_CBP)/sqrt(N) flip(mean(PSTH_rec_CBP)-std(PSTH_rec_CBP)/sqrt(N))],[0 0 1],'FaceAlpha',0.2,'EdgeColor','none');
fill([time flip(time)],[mean(PSTH_sim_CBP)+std(PSTH_sim_CBP)/sqrt(N) flip(mean(PSTH_sim_CBP)-std(PSTH_sim_CBP)/sqrt(N))],[0 0 1],'FaceAlpha',0.5,'EdgeColor','none');
c=plot(time,mean(PSTH_rec),'--','Color',[1 0 0]);
d=plot(time,mean(PSTH_sim),'LineWidth',2,'Color',[1 0 0]);
plot(time,mean(PSTH_rec_CBP),'--','Color',[0 0 1]);
plot(time,mean(PSTH_sim_CBP),'LineWidth',2,'Color',[0 0 1]);
legend([a b c d],'Vision spikes','CBP spikes','Recorded','Simulated')
xlim([time(1) time(end)])
hold off 

disp(fittedGLM.xvalperformance.glm_normedbits)
disp(fittedGLM_CBP.xvalperformance.glm_normedbits)

%% FIGURE 3: correlations and such
rec_diff=mean(PSTH_rec)-mean(PSTH_rec_CBP);
sim_diff=mean(PSTH_sim)-mean(PSTH_sim_CBP);
correlation=xcorr(abs(rec_diff)', abs(sim_diff)', 'unbiased');

figure(3)
subplot(2,1,1)
title('Difference between CBP and Vision')
hold on
plot(time,rec_diff,'--','Color',[1 0 0]);
plot(time,sim_diff,'LineWidth',2,'Color',[0 0 1]);
legend('Recorded','Simulated')
subplot(2,1,2)
hold on
title('Cross Correlation')
plot(correlation)

%% FIGURE 4: correlation with FEM?

load('/Volumes/Data/stimuli/movies/eye-movement/AKH_code_and_images/Krauzlisdata.mat')
index=1;
tdrift=zeros(960,2);
time     = dt*[1:bins];

for j=1:8
    temp=interp1([1:1100],epos_h(mod(j-1,110)+1,1:1100),[1:1000/120:1000],'nearest')';
    temp(:,2)=interp1([1:1100],epos_v(mod(j-1,110)+1,1:1100),[1:1000/120:1000],'nearest');
    temp(:,1)=temp(:,1)-temp(1,1);
    temp(:,2)=temp(:,2)-temp(1,2);
    
    tdrift(index:(index+119),:)=round(temp(1:120,:)*200);
    index=index+120;
end

eye_position=sqrt(tdrift(:,1).^2+tdrift(:,2).^2);
tdrift2=diff(tdrift);
eye_movement=sqrt(tdrift2(:,1).^2+tdrift2(:,2).^2);
CBP_diff=mean(PSTH_rec_CBP)-mean(PSTH_sim_CBP);

figure(4)
hold on
plot(2/120:1/120:8,eye_movement/6,'--','Color',[1 0 0]);
plot(time,abs(CBP_diff),'LineWidth',2,'Color',[0 0 1]);
legend('Eye Movement','Simulated vs. Recorded')

%%
figure(1)
hold on
title('PSTH')
fill([time flip(time)],[mean(PSTH_rec)+std(PSTH_rec)/sqrt(N) flip(mean(PSTH_rec)-std(PSTH_rec)/sqrt(N))],[0 0 1],'FaceAlpha',0.4,'EdgeColor','none');
fill([time flip(time)],[mean(PSTH_sim)+std(PSTH_sim)/sqrt(N) flip(mean(PSTH_sim)-std(PSTH_sim)/sqrt(N))],[1 0 0],'FaceAlpha',0.4,'EdgeColor','none');
a=plot(time,mean(PSTH_rec),'LineWidth',2,'Color',[0 0 1]);
b=plot(time,mean(PSTH_sim),'LineWidth',2,'Color',[1 0 0]);
legend([a b], 'Recorded', 'Simulated')
