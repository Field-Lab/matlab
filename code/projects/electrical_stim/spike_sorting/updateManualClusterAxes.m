function updateManualClusterAxes(H, Data)

% refreshes the axes in the gui manualCluster2D.m

% written by Lauren Hruby, SNL-E
% 2008-10-24

% left axes (included traces)
axes(H.axesLeft); cla
set(H.axesLeft,'YLim',[Data.yMin Data.yMax]);

hold on
traceColor = cool(length(Data.includedIndeces));
for i = 1:length(Data.includedIndeces) %#ok<FXUP>
	current = plot(Data.includedTraces(i,:));
	set(findobj(current,'Type','line'), 'Color', traceColor(i,:))
end 
hold off

%right axes (discluded traces)
axes(H.axesRight); cla
set(H.axesRight,'YLim',[Data.yMin Data.yMax]);

hold on
traceColor = cool(length(Data.discludedIndeces));
for i = 1:length(Data.discludedIndeces) %#ok<FXUP>
	current = plot(Data.discludedTraces(i,:));
	set(findobj(current,'Type','line'), 'Color', traceColor(i,:))
end
hold off