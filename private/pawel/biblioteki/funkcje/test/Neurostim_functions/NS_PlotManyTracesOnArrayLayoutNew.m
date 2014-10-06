function h=NS_PlotManyTracesOnArrayLayoutNew(Waveforms,Channels,WaveformTypes,ArrayID,FigureNumber,FigureProperties,NS_GlobalConstants);
%This function does NOT plot traces from file.
%Plots many traces showing raw data for given channel. For ewach trace the
%starting timing point is defined in 'Timings' input array.
%If ArtifatCancellation=1, then if the WaveformTypes(i)=0, it is assumed to
%be artifact (on all channels at the same time). Then artifacts are
%averaged for each channel and subtracted from each trace.

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

TimeRange=FigureProperties.TimeRange;
AmplitudeRange=FigureProperties.AmplitudeRange;
FontSize=FigureProperties.FontSize;
Colors=FigureProperties.Colors;

%t=[TimeRange(1):TimeRange(2)]/Fs*1000;

Xstep=60;
Ystep=60;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

t=[TimeRange(1):TimeRange(2)]; %/Fs*1000;

Coordinates=zeros(2,length(Channels));
ChannelsNumber=length(Channels);
for i=1:ChannelsNumber
    Coordinates(1,i)=electrodeMap.getXPosition(Channels(i));
    Coordinates(2,i)=electrodeMap.getYPosition(Channels(i));
end

X0=min(Coordinates(1,:));
X1=max(Coordinates(1,:));
Y0=min(Coordinates(2,:));
Y1=max(Coordinates(2,:));

a=find(Coordinates(2,:)==Y0);
z=10000;
for i=1:length(a)
    if Coordinates(1,a(i))<z
        z=Coordinates(1,a(i));
        SE=a(i);
    end
end

Xspread=X1-X0+Xstep*1.5;
Yspread=Y1-Y0+Ystep*1.5;
data=Waveforms;
figure(FigureNumber);
%clf;
for i=1:ChannelsNumber    
    X=(Coordinates(1,i)-X0)/Xspread+0.05;
    Y=(Coordinates(2,i)-Y0)/Yspread+0.1;
    XSize=Xstep/Xspread;
    YSize=Ystep/Yspread;
    subplot('Position',[X Y XSize*0.9 YSize*0.9]);
    hold on;
    for j=1:length(Timings)        
        signal=double(data(:,Channels(i)+1)); %first channel is a TTL channel
        %signal=Waveforms;
            if WaveformTypes(j)~=0 %if this is NOT an artifact...
                WaveformColorIndex=WaveformTypes(j)-floor(WaveformTypes(j)/length(Colors))*length(Colors)+1;  
                %signal=signal-Artifact(:,i);
                plot(t,signal,Colors(WaveformColorIndex));                    
                hold on;
            end        
    end            
    %h=text(TimeRange(1)/Fs*1000+(TimeRange(2)/Fs*1000-TimeRange(1)/Fs*1000)*0.9,AmplitudeRange(1)+(AmplitudeRange(2)-AmplitudeRange(1))*0.9,num2str(str(i)));
    %set(h,'FontSize',FontSize);
    xlim(TimeRange/Fs*1000);
    ylim(AmplitudeRange);
    grid on;
    h=gca;
    set(h,'Box','on');
    if i==SE        
        set(h,'FontSize',FontSize);
        xlabel('time [ms]');
        ylabel('output signal [mV]');
    else
        set(h,'XTickLabel',[]);
        set(h,'YTickLabel',[]);
    end
end

y=Coordinates;