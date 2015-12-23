function show_all_cones

global datarun cones myCells
persistent myPlot


if ishandle(myPlot)
    delete(myPlot);
end


myPlot=subplot('position',[0.05 0.2 0.3 0.3]);
set(gca,'DataAspectRatio',[1 1 1])
hold on

plot_rf_summaries(datarun, myCells, 'clear', false, 'label', false,'plot_fits', true, 'fit_color', 'k')

for i=1:length(cones)    
    plot(cones{i}(:,1),cones{i}(:,2),'*','color','r','markersize',5)
end

