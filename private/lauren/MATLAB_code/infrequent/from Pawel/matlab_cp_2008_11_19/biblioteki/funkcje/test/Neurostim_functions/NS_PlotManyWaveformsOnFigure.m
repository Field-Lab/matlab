function h=NS_PlotManyWaveformsOnFigure(Traces,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants);
% NIE SKONCZONE !!!!! 2008/03/09
% Traces - only 2-dimensional here: N*M, where N - number of traces, M -
% length of single trace
% h - pointer to the axes
STraces=size(Traces);

Fs=NS_GlobalConstants.SamplingFrequency;

FigureNumber=FigureProperties.FigureNumber;
TimeRange=FigureProperties.TimeRange;
AmplitudeRange=FigureProperties.AmplitudeRange;
FontSize=FigureProperties.FontSize;
Colors=FigureProperties.Colors;
%Subplot=FigureProperties.Subplot;

%TimeRange=TimeRange./Fs*1000;
t=([0:STraces(2)-1]+TimeRange(1))/Fs*1000;
t=[TimeRange(1):TimeRange(2)]/Fs*1000;

figure(FigureNumber);
%if Subplot(1)~=0
%    subplot(Subplot(1),Subplot(2),Subplot(3));
%end
%clf;

for i=1:STraces(1) %for each signature...                              
    signal0=Traces(i,:);          
    WaveformColorIndex=WaveformTypes(i)-floor(WaveformTypes(i)/length(Colors))*length(Colors)+1;    
    h=plot(t,signal0,Colors(WaveformColorIndex)); 
    set(h,'LineWidth',2);
    hold on;    
    xlim(TimeRange);
    %ylim(AmplitudeRange);
    grid on;    
    h=gca;
    set(h,'FontSize',FontSize);
    xlabel('time [ms]');
    ylabel('output signal [mV]');      
    %set(h,'XTickLabel',[]);
    %set(h,'YTickLabel',[]);            
end