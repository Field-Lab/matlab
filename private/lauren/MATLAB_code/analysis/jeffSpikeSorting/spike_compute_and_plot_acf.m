function [] = spike_compute_and_plot_acf(spikes_times,acf_rng,acf_axis,plot_color)

if ~exist('plot_color')
    plot_color = [0 0 1];
    %disp('no color')
end

%sptr=zeros(max(spikes_times/5),1);sptr(round(spikes_times/5))=1;acf_subset=xcorr(sptr,acf_rng);
%figure(acf_fig);plot(acf_subset(acf_rng+2:acf_rng*2+1)/4)

%compute ISIs in msec
ISIs=diff(spikes_times/20);
%remove ISIs greater outside the desired range
ISIs=ISIs(find(ISIs<acf_rng));
%plot histogram
axes(acf_axis);
hist(ISIs,100);
set(findobj(gca,'Type','patch'),'FaceColor',plot_color,'EdgeColor',plot_color)

%title(sprintf('%d spikes',length(spikes_times)));

% label with number of spikes
xmax = get(gca,'XLim');
ymax = get(gca,'YLim');
set(gca,'XLim',[0 acf_rng]);
set(gca,'YLim',[0 ymax(2)*1.2]);
%ymax = get(gca,'YLim');
text(0.1*xmax(2),1.08*ymax(2),sprintf('%d spikes',length(spikes_times)));

drawnow
