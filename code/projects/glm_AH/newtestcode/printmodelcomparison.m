function printmodelcomparison(fittedGLM1,fittedGLM2,comparename,printname)

dt = fittedGLM1.t_bin;
tstim = fittedGLM1.bins_per_frame * dt;
frames    = fittedGLM1.linearfilters.Stimulus.frame_shifts;
cellinfo    = fittedGLM1.cellinfo;
GLMType = fittedGLM1.GLMType;

if GLMType.TonicDrive
    MU1 = fittedGLM1.linearfilters.TonicDrive.Filter;
    MU2 = fittedGLM2.linearfilters.TonicDrive.Filter;
end
if GLMType.PostSpikeFilter
    PS1 = fittedGLM1.linearfilters.PostSpike.Filter;
    PS2 = fittedGLM2.linearfilters.PostSpike.Filter;
end

K1        = fittedGLM1.linearfilters.Stimulus.Filter;
K1_time1  = fittedGLM1.linearfilters.Stimulus.time_rk1; 
K1_space1 = fittedGLM1.linearfilters.Stimulus.space_rk1; 

K2        = fittedGLM2.linearfilters.Stimulus.Filter;
K2_time1  = fittedGLM2.linearfilters.Stimulus.time_rk1; 
K2_space1 = fittedGLM2.linearfilters.Stimulus.space_rk1; 

xval1 = fittedGLM1.xvalperformance.glm_normedbits;
xval2 = fittedGLM2.xvalperformance.glm_normedbits;

obj1 = abs(fittedGLM1.rawfit.objective_val);
obj2 = abs(fittedGLM2.rawfit.objective_val);

%%
clf


subplot(9,1,1)
axis off
set(gca, 'fontsize', 12)
c = 0;
text(0, 2-0.1*c,sprintf('%s: %s %d: %s',cellinfo.exp_nm, cellinfo.celltype,cellinfo.cid, comparename),'interpreter','none')
c = c+4;
text(0, 2-0.1*c,sprintf('In Red: %s',fittedGLM1.GLMType.fitname),'interpreter','none')
c = c+4;
text(0, 2-0.1*c,sprintf('In Blue: %s',fittedGLM2.GLMType.fitname),'interpreter','none')
c = c+4;
text(0, 2-0.1*c,sprintf('(Blue-Red)/Red in Xval score: %1.2e ',(xval2-xval1)/xval1 ),'interpreter','none')
c = c+4;
text(0, 2-0.1*c,sprintf('(Blue-Red)/Red in Maximized Objective Function: %1.2e',(obj2-obj1)/obj1),'interpreter','none')
if isfield(fittedGLM2, 'pt_nonlinearity_param') && ~isfield(fittedGLM1,'pt_nonlinearity_param')
    c = c+4;
    text(0, 2-0.1*c,sprintf('Nonlinearity Param for Blue is %1.2e',fittedGLM2.pt_nonlinearity_param),'interpreter','none')
end




% Space Filter 1
if ~GLMType.CouplingFilters, subplot(9,3,[4 7]); end  % label in 50 msec intervals
set(gca, 'fontsize', 10);
imagesc(K1_space1);
title('Space Filters'); axis off

% Space Filter 2
if ~GLMType.CouplingFilters, subplot(9,3,[10 13]); end  % label in 50 msec intervals
set(gca, 'fontsize', 8);
imagesc(K2_space1);
axis off

 
% Time Filter 1 
LW = 2;
if ~GLMType.CouplingFilters,subplot(9,3,[5 8]) ; end
set(gca, 'fontsize', 10);
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K1_time1,'color', 'r','linewidth', LW); 
xlim([0, time_msec(end)])
title('Time Filters');
set(gca, 'xtick',[])
set(gca, 'ytick', [0]); 

% Time Filter 2
LW = 2;
if ~GLMType.CouplingFilters,subplot(9,3,[11 14]) ; end
set(gca, 'fontsize', 10);
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K2_time1,'color', 'b','linewidth', LW); 
xlim([0, time_msec(end)])
xlabel('msec'); 
set(gca, 'ytick', [0]); 


% Post Spike Filter 1
if ~GLMType.CouplingFilters, subplot(9,3,[6,9]); end
set(gca, 'fontsize', 10);
bins    = [1:length(PS1)];
time_msec = 1000*dt*bins ;
oneline = ones(1,length(time_msec)); 
plot(time_msec, oneline, 'k-'); hold on
plot(time_msec, exp(PS1),'color', 'r','linewidth', LW); 
xlim([0, time_msec(end)]);
ylim([0, max(1.5, max(exp(PS1)))]);
ylabel('gain');  title('Post Spike Filters'); set(gca,'xtick',[])

% Post Spike Filter 2
if ~GLMType.CouplingFilters, subplot(9,3,[12 15]); end
set(gca, 'fontsize', 10);
bins    = [1:length(PS2)];
time_msec = 1000*dt*bins ;
oneline = ones(1,length(time_msec)); 
plot(time_msec, oneline, 'k-'); hold on
plot(time_msec, exp(PS2),'color', 'b','linewidth', LW); 
xlim([0, time_msec(end)]);
ylim([0, max(1.5, max(exp(PS2)))]);
ylabel('gain'); xlabel('msec'); 


%%
subplot(27,1,[19:23])
set(gca, 'fontsize', 10);
bins_per_rate = 10;
secs     = 5;
bins     = 120 * secs * fittedGLM1.bins_per_frame;
trials   = size(fittedGLM1.xvalperformance.rasters.recorded,1);
time     = dt*[1:bins];
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([1 , 3*trials]); hold on
set(gca,'xtick',[0:1:secs]); 
for i_test = [1,2,3,1]
    
    if i_test == 1, xval_rast = fittedGLM1.xvalperformance.rasters.recorded; end
    if i_test == 2, xval_rast = fittedGLM2.xvalperformance.rasters.glm_sim;  end
    if i_test == 3, xval_rast = fittedGLM1.xvalperformance.rasters.glm_sim;  end
    
	rate = sum(xval_rast,1);
    for i_shift = -floor(bins_per_rate/2): ceil(bins_per_rate/2)
        rate = rate + circshift(rate, [0 i_shift]);
    end
    
    rate = rate/bins_per_rate;
    
    if i_test == 1, plot(time,rate,'k'); end
    if i_test == 2, plot(time,rate,'b'); end
    if i_test == 3, plot(time,rate,'r'); end       
end 


%%
subplot(27,1,[24:27])
set(gca, 'fontsize', 10);
secs     = 5;
bins     = 120 * secs * fittedGLM1.bins_per_frame;
trials   = size(fittedGLM1.xvalperformance.rasters.recorded,1);
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([1 , 3*trials]); hold on
set(gca,'xtick',[0:1:secs]); 
set(gca,'ytick',[1, (trials+1), 2*trials + 1]); set(gca, 'yticklabel',{'EXP','GLM2','GLM1'});
xlabel('seconds'); 
for i_test = [1,2,3,1]
    
    if i_test == 1, xval_rast = fittedGLM1.xvalperformance.rasters.recorded; end
    if i_test == 2, xval_rast = fittedGLM2.xvalperformance.rasters.glm_sim;  end
    if i_test == 3, xval_rast = fittedGLM1.xvalperformance.rasters.glm_sim;  end
    
	
    offset = (i_test-1)*trials;
    
    for i_trial = 1:trials
        spikebin = find(xval_rast(i_trial,:));
        if length(spikebin) < (4/trials)*length(find(fittedGLM1.xvalperformance.rasters.recorded));         
            if i_test == 1, plot(dt*spikebin, i_trial*ones(1,length(spikebin)) + offset, 'k.'); end
            if i_test == 2, plot(dt*spikebin, i_trial*ones(1,length(spikebin)) + offset, 'b.'); end
            if i_test == 3, plot(dt*spikebin, i_trial*ones(1,length(spikebin)) + offset, 'r.'); end
        end
    end    
end 
orient landscape
eval(sprintf('print -dpdf %s.pdf',printname))




end