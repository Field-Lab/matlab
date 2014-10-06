function spike_plot_average_spike(spikes,num_electrodes,points_per_electrode,avg_spike_axis,plot_color)

if ~exist('plot_color')
    plot_color = 'b';
    disp('no color')
end

% compute mean
mean_spike = mean(spikes,1);

% use provided axes
axes(avg_spike_axis);

% plot waveform
plot([1 length(mean_spike)],[0 0],'k')
hold on
for ee = 1:num_electrodes
    min_ = (ee-1)*points_per_electrode+1;
    max_ = min_ + points_per_electrode - 1;
    plot(min_:max_,mean_spike(min_:max_),'color',plot_color,'LineWidth',1);
end;

for ee = 1:(num_electrodes-1)
    plot((points_per_electrode*ee+0.5)*[1 1],get(gca,'YLim'),'k')
end

hold off

%xlabel('average spike');
drawnow
