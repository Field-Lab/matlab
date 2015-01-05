function plotraster(xval, fittedGLM, labels, separate, PSTH)

dt = fittedGLM.t_bin;
bins     = 120 * 48 * fittedGLM.bins_per_frame;
time     = dt*[1:bins];
rec_rast = xval.rasters.recorded(:,1:bins);
sim_rast = xval.rasters.glm_sim(:,1:bins);
trials   = size(rec_rast,1);

if separate
    hFig1=figure;
    set(hFig1, 'Position', [100 100 800 250])
    hold on
    for i_trial = 1:trials
        rec1 = time(find(rec_rast(i_trial,:)));
        plot(rec1, i_trial, 'k.')
    end
    ylim([1 trials])
    if labels
        xlabel('Time (seconds)')
        ylabel('Trials')
    else
        axis off
    end
    
    hFig2=figure;
    set(hFig2, 'Position', [100 100 800 250])
    hold on
    for i_trial = 1:trials
        sim1 = time(find(sim_rast(i_trial,:)));
        if length(sim1) < 4*length(rec1)
            plot(sim1, i_trial, 'r.')
        end
    end
    ylim([1 trials])
    if labels
        xlabel('Time (seconds)')
        ylabel('Trials')
    else
        axis off
    end
    
else
    

    hFig=figure;
    
    if PSTH
        set(hFig, 'Position', [100 100 1000 500])
        subplot(2,1,1)
    else
        set(hFig, 'Position', [100 100 1000 250])
    end
    
    hold on
    runawaytrial=ones(trials,1);
    for i_trial = 1:trials
        rec1 = time(find(rec_rast(i_trial,:)));
        plot(rec1, i_trial, 'k.')
        sim1 = time(find(sim_rast(i_trial,:)));
        if length(sim1) < 4*length(rec1)
            plot(sim1, i_trial + trials, 'r.')
                else
        runawaytrial(i_trial)=0;
        end
    end
    ylim([1 trials*2])
    xlim([time(1) time(end)])
    
    if labels
        xlabel('Time (seconds)')
        set(gca,'YTickLabel',{'1','Data',num2str(trials),'Model'},'YTick',[1 trials/2 trials 3*trials/2])
        ylabel('Trials')
    end
   
    if PSTH
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
        
        subplot(2,1,2)
        hold on
        % fill([time flip(time)],[mean(PSTH_rec)+std(PSTH_rec)/sqrt(N) flip(mean(PSTH_rec)-std(PSTH_rec)/sqrt(N))],'k','FaceAlpha',0.3,'EdgeColor','none');
        % fill([time flip(time)],[mean(PSTH_sim)+std(PSTH_sim)/sqrt(N) flip(mean(PSTH_sim)-std(PSTH_sim)/sqrt(N))],'r','FaceAlpha',0.3,'EdgeColor','none');
        plot(time,mean(PSTH_rec),'k');
        plot(time,mean(PSTH_sim),'r');
        xlim([time(1) time(end)])
    end
    
end



end