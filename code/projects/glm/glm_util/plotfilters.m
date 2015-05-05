function plotfilters(fittedGLM)

% Plots the GLM filters

GLMType = fittedGLM.GLMType;
dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;
if GLMType.TonicDrive
    MU = fittedGLM.linearfilters.TonicDrive.Filter;
end
if GLMType.PostSpikeFilter
    PS = fittedGLM.linearfilters.PostSpike.Filter;
end
% NBCoupling
if GLMType.CouplingFilters
    CP = fittedGLM.linearfilters.Coupling.Filter;
end
% end NBCoupling

if strcmp(fittedGLM.GLMType.stimfilter_mode, 'rk2')
    rows = 2;
    pos = [100 100 800 400];
    K_time2  = fittedGLM.linearfilters.Stimulus.time_rk2; 
    K_space2 = fittedGLM.linearfilters.Stimulus.space_rk2; 
else
    rows = 1;
    pos = [100 100 800 200];
end

K        = fittedGLM.linearfilters.Stimulus.Filter;
K_time1  = fittedGLM.linearfilters.Stimulus.time_rk1; 
K_space1 = fittedGLM.linearfilters.Stimulus.space_rk1; 

subplot(rows,4,1)  % label in 50 msec intervals
set(gca, 'fontsize', 10);
imagesc(K_space1);
colormap gray
title('Space Filter'); axis off

LW = 2;
subplot(rows,4,2)
set(gca, 'fontsize', 10);
frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K_time1,'color', 'b','linewidth', LW); 
xlim([time_msec(1), time_msec(end)])
xlabel('msec'); title('Time Filter');
set(gca, 'ytick', [0]); 

if rows == 2
    subplot(rows,4,5)  % label in 50 msec intervals
    set(gca, 'fontsize', 10);
    imagesc(K_space2);
    colormap gray
    title('Space Filter'); axis off
    
    LW = 2;
    subplot(rows,4,6)
    set(gca, 'fontsize', 10);
    frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
    time_msec = 1000*tstim*frames ;
    zeroline = zeros(1,length(time_msec));
    plot(time_msec, zeroline, 'k-'); hold on
    plot(time_msec,K_time2,'color', 'b','linewidth', LW);
    xlim([time_msec(1), time_msec(end)])
    xlabel('msec'); title('Time Filter');
    set(gca, 'ytick', [0]);
end

if GLMType.PostSpikeFilter
    subplot(rows,4,3)
    set(gca, 'fontsize', 10);
    bins    = [1:length(PS)];
    time_msec = 1000*dt*bins ;
    oneline = ones(1,length(time_msec)); 
    plot(time_msec, oneline, 'k-'); hold on
    plot(time_msec, exp(PS),'color', 'b','linewidth', LW);
    xlim([0, time_msec(end)]);
    ylim([0, max(1.5, max(exp(PS)))]);
    ylabel('gain'); xlabel('msec'); title('Post Spike Filter')
end


% NBCoupling
if GLMType.CouplingFilters
    subplot(rows,4,4)
    set(gca, 'fontsize', 10);
    plot(time_msec, oneline, 'k-'); hold on
    bins    = [1:length(CP{1})];
    time_msec = 1000*dt*bins ;
    for pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
        plot(time_msec, exp(CP{pair}));
    end
    xlim([0, time_msec(end)]);
    ylim([0, 2]);
    ylabel('gain'); xlabel('msec'); title('Coupling Filters')
end
% end NBCoupling
set(gcf, 'Position', pos)

end

    



     