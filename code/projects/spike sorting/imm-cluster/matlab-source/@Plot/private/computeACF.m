function [] = computeACF(spikes_times,acf_rng,acf_axis,plot_color)

if nargin < 4
    plot_color = [0 0 1];
end

%compute ISIs in msec
ISIs=diff(spikes_times/20);
%remove ISIs greater outside the desired range
ISIs=ISIs(find(ISIs<acf_rng));
%plot histogram
axes(acf_axis);
hist(ISIs,100);
set(findobj(gca,'Type','patch'),'FaceColor',plot_color,'EdgeColor',plot_color)

% label with number of spikes
xmax = get(gca,'XLim');
ymax = get(gca,'YLim');
set(gca,'XLim',[0 acf_rng]);
set(gca,'YLim',[0 ymax(2)*1.2]);
%ymax = get(gca,'YLim');
text(0.1*xmax(2),1.08*ymax(2),sprintf('%d spikes',length(spikes_times)));

drawnow