clear

load('/Volumes/Analysis/nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_noCP_CBP_p8IDp8/standardparams/NSEM_mapPRJ/2012-08-09-3/ONPar_1772.mat')
fittedGLM_CBP=fittedGLM;
load('/Volumes/Analysis/nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/fit2/standardparams/NSEM_mapPRJ/2012-08-09-3/ONPar_1772.mat')

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

runawaytrial=zeros(trials,1);
for i_trial = 1:trials
    sim1{i_trial} = time(find(sim_rast(i_trial,:)));
    rec1{i_trial} = time(find(rec_rast(i_trial,:)));
    plot(rec1{i_trial}, i_trial, 'k.')
    if length(sim1{i_trial}) < 4*length(rec1{i_trial}) 
        plot(sim1{i_trial}, i_trial + trials, 'r.')
        runawaytrial(i_trial)=1;
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
    sim2{i_trial} = time(find(sim_rast(i_trial,:)));
    rec2{i_trial} = time(find(rec_rast(i_trial,:)));
    plot(rec2{i_trial}, i_trial, 'k.')
    if length(sim2{i_trial}) < 4*length(rec2{i_trial}) 
        plot(sim2{i_trial}, i_trial + trials, 'r.')
        runawaytrial(i_trial)=1;
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
    i_trial=runtrials(i);
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
mean(mean(PSTH_rec(:,2000:leftover)));
mean(mean(PSTH_sim(:,2000:leftover)));
mean(mean(PSTH_rec_CBP(:,2000:leftover)));
mean(mean(PSTH_sim_CBP(:,2000:leftover)));
hold off

disp('Distance Time')
%% victor distance
cost=5;
cluster_dist_sim=zeros(trials);
for i=1:trials
    spk1=sim1{i};
    for j=1:i
        spk2=sim1{j};
        if length(spk1)<200 && length(spk2)<200
            cluster_dist_sim(i,j)=victor_dist(spk1,spk2,cost);
        else
            cluster_dist_sim(i,j)=0;
        end
    end
    disp(i)
end

disp('1/6')
    
cluster_dist_rec=zeros(trials);
for i=1:trials
    spk1=rec1{i};
    for j=1:i
        spk2=rec1{j};
        if length(spk1)<200 && length(spk2)<200
            cluster_dist_rec(i,j)=victor_dist(spk1,spk2,cost);
        else
            cluster_dist_rec(i,j)=0;
        end
    end
    disp(i)
end


disp('2/6')
    
CBP_dist_sim=zeros(trials);
for i=1:trials
    spk1=sim2{i};
    for j=1:i
        spk2=sim2{j};
        if length(spk1)<200 && length(spk2)<200
            CBP_dist_sim(i,j)=victor_dist(spk1,spk2,cost);
        else
            CBP_dist_sim(i,j)=0;
        end
    end
    disp(i)
end

disp('3/6')
    
CBP_dist_rec=zeros(trials);
for i=1:trials
    spk1=rec2{i};
    for j=1:i
        spk2=rec2{j};
        if length(spk1)<200 && length(spk2)<200
            CBP_dist_rec(i,j)=victor_dist(spk1,spk2,cost);
        else
            CBP_dist_rec(i,j)=0;
        end
    end
    disp(i)
end

disp('4/6')

cross_dist_rec=zeros(trials);
for i=1:trials
    spk1=rec1{i};
    for j=1:trials
        spk2=rec2{j};
        if length(spk1)<200 && length(spk2)<200
            cross_dist_rec(i,j)=victor_dist(spk1,spk2,cost);
        else
            cross_dist_rec(i,j)=0;
        end
    end
    disp(i)
end

disp('5/6')

cross_dist_sim=zeros(trials);
for i=1:trials
    spk1=sim1{i};
    for j=1:trials
        spk2=sim2{j};
        if length(spk1)<200 && length(spk2)<200
            cross_dist_sim(i,j)=victor_dist(spk1,spk2,cost);
        else
            cross_dist_sim(i,j)=0;
        end
    end
    disp(i)
end
    
    
    
    
    
    
    
    
    
    
    



