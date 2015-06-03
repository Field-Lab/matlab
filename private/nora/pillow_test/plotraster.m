function plotraster(xval, fittedGLM, varargin)
% Optional Parameters
p=inputParser;
addRequired(p,'xval',@isstruct);
addRequired(p,'fittedGLM', @isstruct);
addParameter(p,'raster_length',5, @isnumeric); % length of raster in seconds, default 5
addParameter(p,'labels',false, @islogical); % display axis ticks and labels if true
addParameter(p,'separate',false, @islogical); % separate windows (true) or in subplots (false)
addParameter(p,'PSTH',true, @islogical); % display PSTH if true
addParameter(p,'start_time',0,@isnumeric); % set the raster start time to something other than 0
parse(p,xval, fittedGLM,varargin{:});

if ~isfield(xval.rasters,'recorded')
    predict_only=true;
    rec1=zeros(10^6,1);
else
    predict_only=false;
end

dt = fittedGLM.t_bin;
bins     = fittedGLM.GLMPars.approx_refresh_hz * p.Results.raster_length * fittedGLM.bins_per_frame;
start_bin = fittedGLM.GLMPars.approx_refresh_hz * p.Results.start_time * fittedGLM.bins_per_frame + 1;
bins_idx=start_bin:(start_bin+bins-1);
time     = dt*bins_idx;
try
    if ~predict_only
        rec_rast = xval.rasters.recorded(:,bins_idx);
    end
    sim_rast = xval.rasters.glm_sim(:,bins_idx);
catch
    error('The rasters are not long enough')
end
trials   = size(sim_rast,1);

if p.Results.separate
    
    if ~predict_only
        hFig1=figure;
        set(hFig1, 'Position', [100 100 800 250])
        hold on
        for i_trial = 1:trials
            rec1 = time(find(rec_rast(i_trial,:)));
            plot(rec1, i_trial, 'k.')
        end
        ylim([1 trials])
        if p.Results.labels
            xlabel('Time (seconds)')
            ylabel('Trials')
        else
            axis off
        end
    end
    
    hFig2=figure;
    set(hFig2, 'Position', [100 100 800 250])
    hold on
    runawaytrial=ones(trials,1);
    for i_trial = 1:trials
        sim1 = time(find(sim_rast(i_trial,:)));
        if length(sim1) < 4*length(rec1)
            plot(sim1, i_trial, 'r.')
        else
            runawaytrial(i_trial)=0;
        end
    end
    ylim([1 trials])
    if p.Results.labels
        xlabel('Time (seconds)')
        ylabel('Trials')
    else
        axis off
    end
    
else
    
    hFig=figure;
    
    if p.Results.PSTH
        set(hFig, 'Position', [100 100 1500 500])
        subplot(2,1,1)
    else
        set(hFig, 'Position', [100 100 1500 250])
    end
    
    hold on
    runawaytrial=ones(trials,1);
    for i_trial = 1:trials
        if ~predict_only
            rec1 = time(find(rec_rast(i_trial,:)));
            plot(rec1, i_trial, 'k.')
            sim1 = time(find(sim_rast(i_trial,:)));
            if length(sim1) < 4*length(rec1)
                plot(sim1, i_trial + trials, 'r.')
            else
                runawaytrial(i_trial)=0;
            end
            ylim([1 trials*2])
        else
            sim1 = time(find(sim_rast(i_trial,:)));
            if length(sim1) < 4*length(rec1)
                plot(sim1, i_trial, 'r.')
            else
                runawaytrial(i_trial)=0;
            end
            ylim([1 trials])
        end
    end
    
    xlim([time(1) time(end)])
    
    if p.Results.labels
        xlabel('Time (seconds)')
        set(gca,'YTickLabel',{'1','Data',num2str(trials),'Model'},'YTick',[1 trials/2 trials 3*trials/2])
        ylabel('Trials')
    end
    
end

if p.Results.PSTH
    runtrials=find(runawaytrial);
    N=length(runtrials);
    convolve=30;
    PSTH_rec=zeros(length(runtrials),bins);
    PSTH_sim=zeros(length(runtrials),bins);
    
    for i=1:length(runtrials)
        i_trial=runtrials(i);
        if ~predict_only
            PSTH_rec(i,:)=conv(rec_rast(i_trial,:),ones(convolve,1),'same');
        end
        PSTH_sim(i,:)=conv(sim_rast(i_trial,:),ones(convolve,1),'same');
    end
    
    if p.Results.separate
        hFig3=figure;
        set(hFig3, 'Position', [100 100 800 250])
    else
        subplot(2,1,2)
    end
    hold on
    % fill([time flip(time)],[mean(PSTH_rec)+std(PSTH_rec)/sqrt(N) flip(mean(PSTH_rec)-std(PSTH_rec)/sqrt(N))],'k','FaceAlpha',0.3,'EdgeColor','none');
    % fill([time flip(time)],[mean(PSTH_sim)+std(PSTH_sim)/sqrt(N) flip(mean(PSTH_sim)-std(PSTH_sim)/sqrt(N))],'r','FaceAlpha',0.3,'EdgeColor','none');
    if ~predict_only; plot(time,mean(PSTH_rec),'k'); end
    plot(time,mean(PSTH_sim),'r');
    xlim([time(1) time(end)])
end


end