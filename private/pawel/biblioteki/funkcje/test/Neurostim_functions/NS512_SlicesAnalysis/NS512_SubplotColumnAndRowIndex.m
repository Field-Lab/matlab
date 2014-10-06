function [ColumnIndex,RowIndex]=NS512_SubplotColumnAndRowIndex(NumberOfPlots,PlotIndex);

NumberOfColumns=ceil(sqrt(NumberOfPlots))
NumberOfRows=ceil(NumberOfPlots/NumberOfColumns)

ColumnIndex=floor(PlotIndex/NumberOfColumns)
RowIndex=PlotIndex-(ColumnIndex-1)*NumberOfColumns