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



%{
% Space Filter 1
if ~GLMType.CouplingFilters, subplot(9,3,[4]); end  % label in 50 msec intervals
set(gca, 'fontsize', 10);
imagesc(K1_space1);
title('Space Filters'); axis off

% Space Filter 2
if ~GLMType.CouplingFilters, subplot(9,3,[7]); end  % label in 50 msec intervals
set(gca, 'fontsize', 8);
imagesc(K2_space1);
axis off
%}
 
% Time Filter 1 
LW = 2;
if ~GLMType.CouplingFilters,subplot(9,4,[5 9]) ; end
set(gca, 'fontsize', 10);
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K1_time1,'color', 'r','linewidth', LW); 
xlim([0, time_msec(end)])
title('Time Filter 1');
 xlabel('msec'); 
set(gca, 'ytick', [0]); 

% Time Filter 2
LW = 2;
if ~GLMType.CouplingFilters, subplot(9,4,[7 11]); end
set(gca, 'fontsize', 10);
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K2_time1,'color', 'b','linewidth', LW); 
xlim([0, time_msec(end)])
xlabel('msec'); title('Time Filter 2');
set(gca, 'ytick', [0]); 


% Post Spike Filter 1
if ~GLMType.CouplingFilters,subplot(9,4,[6 10]) ; end
set(gca, 'fontsize', 10);
bins    = [1:length(PS1)];
time_msec = 1000*dt*bins ;
oneline = ones(1,length(time_msec)); 
plot(time_msec, oneline, 'k-'); hold on
plot(time_msec, exp(PS1),'color', 'r','linewidth', LW); 
xlim([0, time_msec(end)]);
ylim([0, max(1.5, max(exp(PS1)))]);
ylabel('gain'); xlabel('msec');   title('Post Spike Filter 1'); %set(gca,'xtick',[])

% Post Spike Filter 2
if ~GLMType.CouplingFilters, subplot(9,4,[8 12]); end
set(gca, 'fontsize', 10);
bins    = [1:length(PS2)];
time_msec = 1000*dt*bins ;
oneline = ones(1,length(time_msec)); 
plot(time_msec, oneline, 'k-'); hold on
plot(time_msec, exp(PS2),'color', 'b','linewidth', LW); 
xlim([0, time_msec(end)]);
ylim([0, max(1.5, max(exp(PS2)))]); title('Post Spike Filter 1');
ylabel('gain'); xlabel('msec'); 


%%%  Raster based stuff
secs = 2;
bins_per_rate = 12;

subplot(10,1,[7:8])
set(gca, 'fontsize', 10);


offset_seconds = .5;
offset_bins    = round(offset_seconds/dt);

bins     = 120 * secs * fittedGLM1.bins_per_frame;
trials   = size(fittedGLM1.xvalperformance.rasters.recorded,1);
time     = dt*[1:bins]; hold on; title(sprintf('Instantaneous Firing Rate with %d msec window', round(1000*bins_per_rate/1200)))
xlim([0, ceil(time(end))]); ylabel('hertz')
set(gca,'xtick',[]); 
for i_test = [2,3,1]
    
    if i_test == 1, xval_rast = fittedGLM1.xvalperformance.rasters.recorded; end
    if i_test == 2, xval_rast = fittedGLM2.xvalperformance.rasters.glm_sim;  end
    if i_test == 3, xval_rast = fittedGLM1.xvalperformance.rasters.glm_sim;  end
    
	rate = sum(xval_rast,1);
    rate0 = sum(xval_rast,1);
    for i_shift = -floor(bins_per_rate/2): ceil(bins_per_rate/2)
        rate = rate + circshift(rate0, [0 i_shift]);
    end
    
    rate            = rate/(bins_per_rate*trials);
    rate_hz_10msec  = rate/(dt); 
    if i_test == 1, spikerate = rate_hz_10msec; end
    if i_test == 2, glmrate_2 = rate_hz_10msec; end
    if i_test == 3, glmrate_1 = rate_hz_10msec; end
    if i_test == 1, plot(time,rate_hz_10msec((1+offset_bins):(bins+offset_bins)),'color','k','linewidth',1);ylim([0 , max(rate_hz_10msec)]); end  
    if i_test == 2, plot(time,rate_hz_10msec((1+offset_bins):(bins+offset_bins)),'color','b','linewidth',.3);ylim([0 , max(rate_hz_10msec)]); end
    if i_test == 3, plot(time,rate_hz_10msec((1+offset_bins):(bins+offset_bins)),'color','r','linewidth',.3);ylim([0 , max(rate_hz_10msec)]); end
end 


% ERROR
subplot(20,1,[18:20])
set(gca, 'fontsize', 10);
bins     = 120 * secs * fittedGLM1.bins_per_frame;
trials   = size(fittedGLM1.xvalperformance.rasters.recorded,1);
time     = dt*[1:bins]; hold on;
xlim([0, ceil(time(end))]);
set(gca,'xtick',[0:1:secs]); title('Error in Firing Rate')
xlabel('SECONDS'); ylabel('hertz')
plot(time, zeros(size(time)), 'k');
for i_model = 1:2
    

    if i_model == 1
        rate_diff = -(spikerate - glmrate_1);
        plot(time,rate_diff((1+offset_bins):(bins+offset_bins)),'color','r'); 
    end
    
    if i_model == 2
        rate_diff = -(spikerate - glmrate_2);
        plot(time,rate_diff((1+offset_bins):(bins+offset_bins)),'color','b'); 
    end
end



%{
subplot(9,1,[9])
set(gca, 'fontsize', 10);


bins     = 120 * secs * fittedGLM1.bins_per_frame;
trials   = size(fittedGLM1.xvalperformance.rasters.recorded,1);
time     = dt*[1:bins]; hold on;
xlim([0, ceil(time(end))]);
set(gca,'xtick',[0:1:secs]); 
xlabel('SECONDS')
plot(time, zeros(size(time)), 'k'); hold on
l2_error_dif = ((glmrate_2 - spikerate).^2) - ((glmrate_1 - spikerate).^2);
plot(time, l2_error_dif(1:bins), 'm')
%}



% RASTERS
subplot(18,1,[8:10]); 
set(gca, 'fontsize', 10);
title('Rasters')
bins     = 120 * secs * fittedGLM1.bins_per_frame;
trials   = size(fittedGLM1.xvalperformance.rasters.recorded,1);
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([1 , 3*trials]); hold on
set(gca,'xtick',[0:1:secs]); 
set(gca,'ytick',[1, (trials+1), 2*trials + 1]); set(gca, 'yticklabel',{'EXP','GLM2','GLM1'});
%xlabel('seconds'); 
for i_test = [1,2,3,1]
    
    if i_test == 1, xval_rast = fittedGLM1.xvalperformance.rasters.recorded; end
    if i_test == 2, xval_rast = fittedGLM2.xvalperformance.rasters.glm_sim;  end
    if i_test == 3, xval_rast = fittedGLM1.xvalperformance.rasters.glm_sim;  end
    
	
    offset = (i_test-1)*trials;
    
    for i_trial = 1:trials
        spikebin = find(xval_rast(i_trial,:)) ;%- offset_seconds;
        if length(spikebin) < (4/trials)*length(find(fittedGLM1.xvalperformance.rasters.recorded));         
            if i_test == 1, plot(dt*spikebin-offset_seconds, i_trial*ones(1,length(spikebin)) + offset, 'k.'); end
            if i_test == 2, plot(dt*spikebin-offset_seconds, i_trial*ones(1,length(spikebin)) + offset, 'b.'); end
            if i_test == 3, plot(dt*spikebin-offset_seconds, i_trial*ones(1,length(spikebin)) + offset, 'r.'); end
        end
    end    
end 
orient landscape
eval(sprintf('print -dpdf %s.pdf',printname))




end