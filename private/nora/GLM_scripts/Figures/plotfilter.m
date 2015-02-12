function plot_filter(fittedGLM, filter_type)

% fittedGLM is the fit you want to get the filter from
% filter_type is 'PS' or 'Temp' for now

dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;

switch filter_type
    case 'PS'
        filter=fittedGLM.linearfilters.PostSpike.Filter;
    case 'Temp'
        filter=fittedGLM.linearfilters.Stimulus.time_rk1;
    case 'CP'
        error('Not yey available')
end

figure
hold on
bins    = [1:length(filter)];
time_msec = 1000*dt*bins ;
oneline = ones(1,length(time_msec)); 
plot(time_msec, oneline, 'k-'); hold on
plot(time_msec,exp(filter),'black','LineWidth',1.5)
  xlim([0, time_msec(end)]);
    ylim([0, max(1.5, max(exp(filter)))]);
    ylabel('Gain'); xlabel('Time (ms)');
    
end