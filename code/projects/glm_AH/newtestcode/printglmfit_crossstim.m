function printglmfit_crossstim(fittedGLM,xval_WN, xval_NSEM,printname)

info    = fittedGLM.cellinfo;
GLMType = fittedGLM.GLMType;
dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;
if GLMType.TonicDrive
    MU = fittedGLM.linearfilters.TonicDrive.Filter;
end
if GLMType.PostSpikeFilter
    PS = fittedGLM.linearfilters.PostSpike.Filter;
end

K        = fittedGLM.linearfilters.Stimulus.Filter;
K_time1  = fittedGLM.linearfilters.Stimulus.time_rk1; 
K_space1 = fittedGLM.linearfilters.Stimulus.space_rk1; 
 

%%
clf


subplot(9,1,1)
axis off
set(gca, 'fontsize', 12)
c = 0;
text(0, 1-0.1*c,sprintf('%s: %s %d: %s-Fit',info.exp_nm, info.celltype,info.cid, GLMType.fit_type))

    

if ~GLMType.CouplingFilters, subplot(9,3,[4 7]); end  % label in 50 msec intervals
set(gca, 'fontsize', 10);
imagesc(K_space1);
title('Space Filter'); axis off

LW = 2;
if ~GLMType.CouplingFilters,subplot(9,3,[13 16]) ; end
set(gca, 'fontsize', 10);
frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K_time1,'color', 'b','linewidth', LW); 
xlim([0, time_msec(end)])
xlabel('msec'); title('Time Filter');
set(gca, 'ytick', [0]); 


if ~GLMType.CouplingFilters, subplot(9,3,[22,25]); end
set(gca, 'fontsize', 10);
bins    = [1:length(PS)];
time_msec = 1000*dt*bins ;
oneline = ones(1,length(time_msec)); 
plot(time_msec, oneline, 'k-'); hold on
plot(time_msec, exp(PS),'color', 'b','linewidth', LW); 
xlim([0, time_msec(end)]);
ylim([0, max(1.5, max(exp(PS)))]);
ylabel('gain'); xlabel('msec'); title('Post Spike Filter')


clear i_test

for i_test = [1,2,1]
    
    if i_test == 1, xval = xval_WN;  teststim = 'WN'; subplot(10,3,[5,6,8,9,11,12,14,15]); end
    if i_test == 2, xval = xval_NSEM; teststim = 'NSEM'; subplot(10,3,(15+[5,6,8,9,11,12,14,15]) ); end
    set(gca, 'fontsize', 10);
	rec_rast = xval.rasters.recorded(:,1:bins);
	sim_rast = xval.rasters.glm_sim(:,1:bins); 
    score    = xval.glm_normedbits;
    
    secs     = 4;
    bins     = 120 * secs * fittedGLM.bins_per_frame;
    trials   = size(rec_rast,1);
    time     = dt*[1:bins];
    
    xlim([0, ceil(time(end))]);
    ylim([1 , 2*trials]); hold on
    for i_trial = 1:trials
        sim1 = time(find(sim_rast(i_trial,:)));
        rec1 = time(find(rec_rast(i_trial,:)));
        plot(rec1, i_trial*ones(1,length(rec1)), 'k.')
        plot(sim1, (i_trial + trials)*ones(1,length(sim1)), 'r.')
    end
    set(gca,'xtick',[0:1:secs]); 
    set(gca,'ytick',[1, (trials+1)]); set(gca, 'yticklabel',{'EXP','GLM'});
    xlabel('seconds'); 
    title(sprintf('Test GLM with %s:  Scored  %d out of 100', teststim, round(100*score) ) );
    
end 

orient landscape
eval(sprintf('print -dpdf %s.pdf', printname))
end