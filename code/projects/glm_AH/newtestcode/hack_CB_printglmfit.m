function hack_printglmfit_CB(fittedGLM,printname)



info    = fittedGLM.cellinfo;
GLMType = fittedGLM.GLMType;
dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;


if isfield(fittedGLM, 'xvalperformance')
    performance = fittedGLM.xvalperformance;
end
if GLMType.TonicDrive
    MU = fittedGLM.linearfilters.TonicDrive.Filter;
end
if GLMType.PostSpikeFilter
    PS = fittedGLM.linearfilters.PostSpike.Filter;
end

K_ex       = fittedGLM.linearfilters.Stimulus.Excitatory_Filter;
K_ex_time  = fittedGLM.linearfilters.Stimulus.Excitatory_Time;
K_ex_space = fittedGLM.linearfilters.Stimulus.Excitatory_Space;


K_in       = fittedGLM.linearfilters.Stimulus.Inhibitory_Filter;
K_in_time  = fittedGLM.linearfilters.Stimulus.Inhibitory_Time;
K_in_space = fittedGLM.linearfilters.Stimulus.Inhibitory_Space;


klen = sqrt(length(K_ex_space(:)));


K_ex_space = reshape(K_ex_space, [klen,klen]);
K_in_space = reshape(K_in_space, [klen,klen]);

clf
subplot(5,3,[1 2])
axis off
set(gca, 'fontsize', 12)
c = 0;
text(-.1, 1-0.1*c,sprintf('%s: %s %d',info.exp_nm, info.celltype,info.cid))
c = c + 1.5;
text(-.1, 1-0.1*c, sprintf('%s-Fit',GLMType.fit_type))
if isfield(fittedGLM, 'xvalperformance')
    c = c + 1.5;
    text(-.1, 1-0.1*c,sprintf('Cross Validated: Bits Per spike %1.3e',performance.logprob_glm_bpspike))
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

LW = 2;
% Plot PS Filter
subplot(5,3,3)
set(gca, 'fontsize', 10);
bins    = [1:length(PS)];
time_msec = 1000*dt*bins ;
oneline = ones(1,length(time_msec)); 
plot(time_msec, oneline, 'k-'); hold on
plot(time_msec, exp(PS),'color', 'b','linewidth', LW); 
xlim([0, time_msec(end)]);
ylim([0, max(1.5, max(exp(PS)))]);
ylabel('gain'); xlabel('msec'); title('Post Spike Filter');

    
% Excitatory Filters
subplot(5,4,[5 9])  % label in 50 msec intervals
set(gca, 'fontsize', 10);
imagesc(K_ex_space);
title('Excite Space'); axis off


subplot(5,4,[6,10])
set(gca, 'fontsize', 10);
frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K_ex_time,'color', 'b','linewidth', LW); 
xlim([time_msec(1), time_msec(end)])
xlabel('msec'); title('Excite Time');
set(gca, 'ytick', [0]); 

% Inhibitory Filters
subplot(5,4,[7,11])  % label in 50 msec intervals
set(gca, 'fontsize', 10);
imagesc(K_in_space);
title('Inhib Space'); axis off

LW = 2;
subplot(5,4,[8,12]);
set(gca, 'fontsize', 10);
frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K_in_time,'color', 'b','linewidth', LW); 
xlim([time_msec(1), time_msec(end)])
xlabel('msec'); title('Inhib Time');
set(gca, 'ytick', [0]); 




%
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
%}    
    
    




orient landscape
eval(sprintf('print -dpdf %s.pdf',printname))
end

    



     