function printglmfit(fittedGLM,printname,optional_note)

info    = fittedGLM.cellinfo;
GLMType = fittedGLM.GLMType;
dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;
performance = fittedGLM.xvalperformance;
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

K        = fittedGLM.linearfilters.Stimulus.Filter;
K_time1  = fittedGLM.linearfilters.Stimulus.time_rk1; 
K_space1 = fittedGLM.linearfilters.Stimulus.space_rk1; 
    
clf
subplot(5,1,1)
axis off
set(gca, 'fontsize', 12)
c = 0;
text(-.1, 1-0.1*c,sprintf('%s: %s %d: %s-Fit',info.exp_nm, info.celltype,info.cid, GLMType.fit_type))
c = c + 1.5;
text(-.1, 1-0.1*c,sprintf('Cross Validated Performance, in Bits Per spike %1.3e',performance.logprob_glm_bpspike))
% AKH 2015-07-15 added room for an extra note
if exist('optional_note','var')
    c = c + 1.5;
    text(-.1, 1-0.1*c,sprintf('%s',optional_note),'interpreter','none');
end

 c = c + 1.5;
text(-.1, 1-0.1*c,sprintf('Fit Type: %s',GLMType.fitname), 'interpreter','none')
c = c + 1.5;
if GLMType.TonicDrive
    text(-.1, 1-0.1*c,sprintf('Tonic drive with gray stim: %1.2e hz',round(exp(MU)))  )
    c = c + 1.5;
end
text(-.1, 1-0.1*c,sprintf('Optimum fmin: %1.4e',fittedGLM.rawfit.objective_val))
c = c + 1.5;
text(-.1, 1-0.1*c,sprintf('Computated at %s',fittedGLM.fit_time))
    

subplot(5,4,[5 9])  % label in 50 msec intervals
set(gca, 'fontsize', 10);
imagesc(K_space1);
title('Space Filter'); axis off

LW = 2;
subplot(5,4,[6,10])
set(gca, 'fontsize', 10);
frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K_time1,'color', 'b','linewidth', LW); 
xlim([time_msec(1), time_msec(end)])
xlabel('msec'); title('Time Filter');
set(gca, 'ytick', [0]); 



if GLMType.PostSpikeFilter
    subplot(5,4,[7,11])
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
    subplot(5,4,[8,12])
    set(gca, 'fontsize', 10);
    
    bins    = [1:length(CP{1})];
    time_msec = 1000*dt*bins ;
    oneline = ones(1,length(time_msec));
    plot(time_msec, oneline, 'k-'); hold on
    for pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
        plot(time_msec, exp(CP{pair}));
    end
    xlim([0, time_msec(end)]);
    ylim([0, 2]);
    ylabel('gain'); xlabel('msec'); title('Coupling Filters')
end
% end NBCoupling

subplot(5,1,[4,5]);

secs     = 8;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
rec_rast = fittedGLM.xvalperformance.rasters.recorded(:,1:bins);
sim_rast = fittedGLM.xvalperformance.rasters.glm_sim(:,1:bins);
trials   = size(rec_rast,1);
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([1 , 2*trials]); hold on

for i_trial = 1:trials
    sim1 = time(find(sim_rast(i_trial,:)));
    rec1 = time(find(rec_rast(i_trial,:)));
    
    
    
    plot(rec1, i_trial, 'k.')
    
    if length(sim1) < 4*length(rec1)
        if length(sim1) > 0
            plot(sim1, i_trial + trials, 'r.')
        end
    end
end

orient landscape
eval(sprintf('print -dpdf %s.pdf',printname))

end

    



     