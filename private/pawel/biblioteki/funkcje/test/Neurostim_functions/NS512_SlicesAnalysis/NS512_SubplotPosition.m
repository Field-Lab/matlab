function Position=NS512_SubplotPosition(NumberOfPlots,PlotIndex,PositionOfWholePlot);
%this function allows you to define soem specific area of the figure
%(PositionOfWholePlot) and plot a number of subplots within this area. The
%function returns the position of specific suplot, calculating this for specific plot based
%on number of plots.

NumberOfColumns=ceil(sqrt(NumberOfPlots))
NumberOfRows=ceil(NumberOfPlots/NumberOfColumns)

ColumnIndex=floor(PlotIndex/NumberOfColumns)
RowIndex=PlotIndex-(ColumnIndex-1)*NumberOfColumns

Position=2;