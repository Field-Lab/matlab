function h=NS_PlotManyTracesOnFigure(FileName,Timings,WaveformTypes,Channel,OffsetCancellation,OffsetSamples,FigureProperties,Delay,NS_GlobalConstants)
%Plots many traces showing raw data for given channel. For ewach trace the
%starting timing point is defined in 'Timings' input array.
%WaveformTypes - there is a specific trace corresponding to each value in
%'Timings' table, but the traces can be of several different types (for
%example, they may correspond to spikes from different neurons). The
%WaveformTypes array, whichc must have the same length as Timings array,
%defines the waveform type for each event. Each waveform type will be
%finally ploted in different color. By default, type 0 corresponds to "no
%spikes" waveform (which is equal to "artifact only"). 

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

FigureNumber=FigureProperties.FigureNumber;
TimeRange=FigureProperties.TimeRange;
AmplitudeRange=FigureProperties.AmplitudeRange;
FontSize=FigureProperties.FontSize;
FontName=FigureProperties.FontName;
Colors=FigureProperties.Colors;

t=[TimeRange(1):TimeRange(2)]/Fs*1000+0.005;
%t=t-0.35;

full_path=[pwd '\' 'data' FileName];
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

figure(FigureNumber);
%subplot(2,1,2);
clf;
if OffsetCancellation==3
    Offset=mean(OffsetSamples);
end

for i=1:length(Timings)
    data=rawFile.getData(Timings(i)+TimeRange(1),TimeRange(2)-TimeRange(1)+1);
    signal=data(:,Channel+1); %first channel is a TTL channel    
    
    switch OffsetCancellation
        case 1
            signal=signal-mean(signal(:,1));
        case 2
            signal=signal-mean(signal(OffsetSamples,1));  
        case 3
            signal=signal-Offset;
    end           
    WaveformColorIndex=WaveformTypes(i)-floor(WaveformTypes(i)/length(Colors))*length(Colors)+1;    
    plot(t,signal/0.84,Colors(WaveformColorIndex)); 
    h=gca;
    set(h,'YLim',AmplitudeRange);
    if Delay==1
        %h=gca;
        set(h,'XLim',[0.6 1.1]);
        %set(h,'YLim',AmplitudeRange);
        i
        refresh;
        pause(2);
    end
    hold on;
end

%axis([min(t) max(t) AmplitudeRange(1) AmplitudeRange(2)]);
h=gca;
%set(h,'YLim',[-331.5 -330.6]);
%set(h,'XLim',[0.734 0.742]);

set(h,'FontSize',FontSize);
set(h,'FontName',FontName);
xlabel('Time [ms]');
ylabel('Signal [\muV]');
grid on;
    
    
