function plotfilters(fittedGLM)

% Plots the GLM filters using the output from glm_fit

GLMType = fittedGLM.GLMType;
dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;
plots = 2; % Number of subplots needed. Start with 2 for space and time

% Get filters
if GLMType.TonicDrive
    MU = fittedGLM.linearfilters.TonicDrive.Filter;
end
if GLMType.PostSpikeFilter
    PS = fittedGLM.linearfilters.PostSpike.Filter;
    plots = plots+1;
    PS_loc = plots;
end
if GLMType.CouplingFilters
    CP = fittedGLM.linearfilters.Coupling.Filter;
    plots = plots+1;
    CP_loc = plots;
end
if GLMType.Subunits
    SU = fittedGLM.SU_filter;
    plots = plots+1;
    SU_loc = plots;
end

if strcmp(fittedGLM.GLMType.stimfilter_mode, 'rk2')
    rows = 2;
    pos = [100 100 200*plots 400];
    K_time2  = fittedGLM.linearfilters.Stimulus.time_rk2; 
    K_space2 = fittedGLM.linearfilters.Stimulus.space_rk2; 
else
    rows = 1;
    pos = [100 100 200*plots 200];
end


% Plot the rank 1 stimulus filters
K        = fittedGLM.linearfilters.Stimulus.Filter;
K_time1  = fittedGLM.linearfilters.Stimulus.time_rk1; 
K_space1 = fittedGLM.linearfilters.Stimulus.space_rk1; 

subplot(rows,plots,1)
set(gca, 'fontsize', 10);
imagesc(K_space1);
colormap gray
title('Space Filter'); axis off

LW = 2;
subplot(rows,plots,2)
set(gca, 'fontsize', 10);
frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K_time1,'linewidth', LW); 
xlim([time_msec(1), time_msec(end)])
xlabel('msec'); title('Time Filter');
set(gca, 'ytick', [0]); 

% Plot the RK2 filters if they exist
if rows == 2
    subplot(rows,plots,5)  % label in 50 msec intervals
    set(gca, 'fontsize', 10);
    imagesc(K_space2);
    colormap gray
    title('Space Filter'); axis off
    
    LW = 2;
    subplot(rows,plots,6)
    set(gca, 'fontsize', 10);
    frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
    time_msec = 1000*tstim*frames ;
    zeroline = zeros(1,length(time_msec));
    plot(time_msec, zeroline, 'k-'); hold on
    plot(time_msec,K_time2,'linewidth', LW);
    xlim([time_msec(1), time_msec(end)])
    xlabel('msec'); title('Time Filter');
    set(gca, 'ytick', [0]);
end

if GLMType.PostSpikeFilter
    subplot(rows,plots,PS_loc)
    set(gca, 'fontsize', 10);
    bins    = [1:length(PS)];
    time_msec = 1000*dt*bins ;
    oneline = ones(1,length(time_msec)); 
    plot(time_msec, oneline, 'k-'); hold on
    plot(time_msec, exp(PS),'linewidth', LW);
    xlim([0, time_msec(end)]);
    ylim([0, max(1.5, max(exp(PS)))]);
    ylabel('gain'); xlabel('msec'); title('Post Spike Filter')
end


% NBCoupling
if GLMType.CouplingFilters
    subplot(rows,plots,CP_loc)
    set(gca, 'fontsize', 10);
    plot(time_msec, oneline, 'k-'); hold on
    bins    = [1:length(CP{1})];
    time_msec = 1000*dt*bins ;
    hold on
    for pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
        plot(time_msec, exp(CP{pair}));
    end
    hold off
    xlim([0, time_msec(end)]);
    ylim([0, 2]);
    ylabel('gain'); xlabel('msec'); title('Coupling Filters')
end

if GLMType.Subunits
    subplot(rows, plots, SU_loc)
    set(gca, 'fontsize', 10);
    imagesc(SU)
    axis image
    axis off
    title('Subunit Filter')

% end NBCoupling
set(gcf, 'Position', pos)

end

    



     