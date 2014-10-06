function y=NS_PlotSignatureOnArrayLayout(Signature,Channels,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants);
%OffsetCancellation - may have three values:
%a) if 0, the trace is plotted directly;
%b) if 1, then the mean value of the whole trace is subtracted before
%plotting;
%c) if 2, then the mean value of all the samples Deined in the OffsetSamples
%array is subtracted before plotting (for examples, if one wants to plot a spike or
%artifact tha stats 20 samples after beginning of each trace, it may be
%reasonable to use the first milisecond to calculate the DC offset).
Fs=NS_GlobalConstants.SamplingFrequency;

FigureNumber=FigureProperties.FigureNumber;
TimeRange=FigureProperties.TimeRange;
AmplitudeRange=FigureProperties.AmplitudeRange;
FontSize=FigureProperties.FontSize;
Colors=FigureProperties.Colors;

Xstep=60; 
Ystep=60;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
TimeRange=TimeRange./Fs*1000;
t=[TimeRange(1):(1/Fs*1000):TimeRange(2)];
SSignature=size(Signature);
t=([1:SSignature(2)]-TimeRange(1))/Fs*1000;

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

figure(FigureNumber);
clf;
for i=1:ChannelsNumber    
    X=(Coordinates(1,i)-X0)/Xspread+0.05;
    Y=(Coordinates(2,i)-Y0)/Yspread+0.06;
    XSize=Xstep/Xspread;
    YSize=Ystep/Yspread;
    subplot('Position',[X Y XSize*0.9 YSize*0.9]);
    signal=Signature(i,:);       
    WaveformColorIndex=WaveformTypes(i)-floor(WaveformTypes(i)/length(Colors))*length(Colors)+1;    
    plot(t,signal,Colors(WaveformColorIndex));                             
    text(TimeRange(1)+(TimeRange(2)-TimeRange(1))*0.8,AmplitudeRange(1)+(AmplitudeRange(2)-AmplitudeRange(1))*0.8,num2str(Channels(i)));
    xlim(TimeRange);
    ylim(AmplitudeRange);
    grid on;    
    h=gca;
    if i==SE        
        set(h,'FontSize',FontSize);
        xlabel('t [ms]');
        ylabel('ADC units');
    else
        set(h,'XTickLabel',[]);
        set(h,'YTickLabel',[]);
    end
end

y=Coordinates;