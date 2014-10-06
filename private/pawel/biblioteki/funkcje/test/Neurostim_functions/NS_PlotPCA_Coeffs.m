function h=NS_PlotPCA_Coeffs(NumberOfClusters,Types,Coeffs,FigureProperties);

FigureNumber=FigureProperties.FigureNumber;
TimeRange=FigureProperties.TimeRange;
AmplitudeRange=FigureProperties.AmplitudeRange;
FontSize=FigureProperties.FontSize;
Colors=FigureProperties.Colors;
Subplot=FigureProperties.Subplot;

if Subplot(1)~=0
    subplot(Subplot(1),Subplot(2),Subplot(3));
end
%clf;
%figure(101);
for i=1:NumberOfClusters
    FindType=find(Types==i);
    WaveformColorIndex=i-floor(i/length(Colors))*length(Colors)+1;
    h=plot(Coeffs(FindType,2),Coeffs(FindType,3),'bd');
    set(h,'Color',Colors(WaveformColorIndex));
    if Colors(WaveformColorIndex)=='g'
        set(h,'Color',[0 0.75 0]);
    end
    set(h,'MarkerSize',10);
    hold on;
end
h=gca;
set(h,'FontSize',FontSize);
grid on;
hold off;
y=i; %for now :)