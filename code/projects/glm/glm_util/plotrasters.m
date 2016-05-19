function [PSTH_sim, PSTH_rec] = plotrasters(xval, fittedGLM, varargin)
% NB 2015-06-18
% DESCRIPTION
% This function takes in xval from glm_predict, fittedGLM from glm_fit, and
% plots rasters and PSTH. Red is predicted and black is recorded.
%
% INPUTS
% 
% REQUIRED
% xval: t, fittehe output structure from glm_predict
% fittedGLM: the output structure from glm_fit
%
% Optional Parameters
%
% raster_length, default 6 seconds. if you plot more than 10 seconds, it
%   gets slow hard to see
% labels, default no labels, if true, you see labels and axis ticks on
%   raster
% separate, default no, separates the rasters and PSTH into separate
%   windows
% PSTH, default no, plots the PSTH
% start_time, default 2, can set to start later in the raster
% PSTH_window_size, default 100, in bins, sets the smoothing of the PSTH
%
p=inputParser;
addRequired(p,'xval',@isstruct);
addRequired(p,'fittedGLM', @isstruct);
addParameter(p,'raster_length',6, @isnumeric); % length of raster in seconds, default 5
addParameter(p,'labels',false, @islogical); % display axis ticks and labels if true
addParameter(p,'separate',false, @islogical); % separate windows (true) or in subplots (false)
addParameter(p,'PSTH',false, @islogical); % display PSTH if true
addParameter(p,'start_time',2,@isnumeric); % set the raster start time to something other than 0
addParameter(p,'PSTH_window_size', 100, @isnumeric); % set the amount smoothing in the PSTH
addParameter(p, 'rasters', 1, @islogical); % set to 0 to only plot the PSTH
parse(p,xval, fittedGLM,varargin{:});


[~,B] = system('hostname');
if strcmp(B(1:6), 'bertha')
    scale = 100;
else
    scale = 100;
end
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

if p.Results.rasters
if p.Results.separate
    
    if ~predict_only
        hFig1=figure;
        set(hFig1, 'Position', [1 1 8 3]*scale)
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
    set(hFig2, 'Position', [1 1 8 3]*scale)
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
        set(hFig, 'Position', [1 1 10 3]*scale)
        subplot(2,1,1)
    else
        set(hFig, 'Position', [1 1 10 3]*scale)
    end
    
    hold on
    runawaytrial=ones(trials,1);
    for i_trial = 1:trials
        if ~predict_only
            rec1 = time(find(rec_rast(i_trial,:)));
            if ~isempty(rec1)
                plot(rec1, i_trial, 'k.')
            end
            sim1 = time(find(sim_rast(i_trial,:)));
if length(sim1) < 4*length(rec1) && ~isempty(sim1)
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
    else
        axis off
    end
    
end
end

if p.Results.PSTH
    convolve=gausswin(p.Results.PSTH_window_size);
    if ~predict_only
        PSTH_rec=conv(sum(rec_rast),convolve,'same');
    end
    PSTH_sim=conv(sum(sim_rast),convolve,'same');
    
    if p.Results.separate || ~p.Results.rasters
        hFig3=figure;
        set(hFig3, 'Position', [1 1 8 3]*scale)
    else
        subplot(2,1,2)
    end
    hold on
    if ~predict_only; plot(time,PSTH_rec,'k'); end
    plot(time,PSTH_sim,'r');
    xlim([time(1) time(end)])
end


end
