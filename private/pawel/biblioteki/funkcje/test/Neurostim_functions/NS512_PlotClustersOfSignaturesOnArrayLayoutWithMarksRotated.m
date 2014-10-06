function y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarksRotated(Traces,Channels,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,MarkedChannels2,t);
% NIE SKONCZONE !!!!! 2008/03/09
%OffsetCancellation - may have three values:
%a) if 0, the trace is plotted directly;
%b) if 1, then the mean value of the whole trace is subtracted before
%plotting;
%c) if 2, then the mean value of all the samples Deined in the OffsetSamples
%array is subtracted before plotting (for examples, if one wants to plot a spike or
%artifact tha stats 20 samples after beginning of each trace, it may be
%reasonable to use the first milisecond to calculate the DC offset).
STraces=size(Traces);
if STraces(2)~=length(Channels)
    error('The size of Traces array does not match length of Channels array');
end

%Channels=find(Channels~=9);
Fs=NS_GlobalConstants.SamplingFrequency;

FigureNumber=FigureProperties.FigureNumber;
TimeRange=FigureProperties.TimeRange;
AmplitudeRange=FigureProperties.AmplitudeRange;
FontSize=FigureProperties.FontSize;
Colors=FigureProperties.Colors;
LineWidth=FigureProperties.LineWidth;
YLabel=FigureProperties.YLabel;

%FontSize=10;

NumberOfTraces=STraces(1);
NumberOfChannels=STraces(2);
TraceLength=STraces(3);

Xstep=60; 
Ystep=60;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
TimeRange=TimeRange./Fs*1000;

%t=[TimeRange(1):(1/Fs*1000):TimeRange(2)];
if nargin < 9
    t=([0:STraces(3)-1])*1000/Fs;
end

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

Xspread=X1-X0+Xstep*1.5; %1.5 !!
Yspread=Y1-Y0+Ystep*1.5;
figure(FigureNumber);
clf;

for c=max(WaveformTypes):-1:min(WaveformTypes)
    TracesToPlot=find(WaveformTypes==c);
    ll=length(TracesToPlot);
    maxWT=max(WaveformTypes);
    minWT=max(min(WaveformTypes),1); % WaveformTyp=0 only for artifact - including artifact does not change color associated with given cluster index
    ColorIndex=c-floor((c-1)/(length(Colors)-1))*(length(Colors)-1)+1;
    if c==0
        ColorIndex=1;
    end
    Color=Colors(ColorIndex);
    disp(['cluster ' num2str(c) ', color ' Color ', ' num2str(ll) ' traces'])
    for i=1:NumberOfChannels                
        Y=(-Coordinates(1,i)+X0)/Xspread*1.05+0.8 %0.15;
        X=(Coordinates(2,i)+Y0)/Yspread/1.05-0.01 %0.1;
        
        %X=(Coordinates(2,i)+Y0)/Yspread/1.1-2.8 %0.1;
        %Y=(-Coordinates(1,i)+X0)/Xspread+0.7 %0.15;

        
        XSize=Xstep/Xspread;
        YSize=Ystep/Yspread;        
        subplot('Position',[X Y XSize*0.57 YSize*1.5]);
        Data=reshape(Traces(TracesToPlot,i,:),length(TracesToPlot),TraceLength);
        h=plot(t,Data');
        set(h,'Color',Color);       
        set(h,'LineWidth',LineWidth);
        %h=text(2.2,150,num2str(Channels(i)));
        grid on;
        h=gca;
        set(h,'YLim',AmplitudeRange);
        set(h,'XLim',TimeRange);
        %disp('glum')
        if Channels(i)==316 %i==SE
            set(h,'FontSize',FontSize);  
            set(h,'XTick',[0:0.5:3]);
            set(h,'XTickLabel',{'' '0' '' '1' '','2',''});
            set(h,'YTick',[-120:60:60]);            
            h=xlabel('Time [ms]');
            set(h,'FontSize',FontSize);            
            %h=ylabel('signal [\muV]');
            h=ylabel(YLabel);
            set(h,'FontSize',FontSize);
        %else if Channels(i)==309
        %        set(h,'XTick',[0:0.5:3]);
        %        set(h,'XTickLabel',{'' '0' '' '1' '','2',''});                
        %        set(h,'YTick',[-120:60:60]);
        %       set(h,'YTickLabel',[]);
        %        h=xlabel('time [ms]');
        %        set(h,'FontSize',FontSize);
            else
            %set(h,'FontSize',FontSize);                
            set(h,'XTickLabel',[]);
            set(h,'XTick',[0:0.5:3]);
            set(h,'YTickLabel',[]);
            set(h,'YTick',[-120:60:60]);
            %set(h,'FontSize',FontSize);
            %h=xlabel(num2str(Channels(i)));
            %h=text(1.7,-70,num2str(Channels(i)));
            %set(h,'FontSize',FontSize);
            %h=xlabel(num2str(Channels(i)));           
            %end
        end   
                                
        hold on;
        %dfg=find(MarkedChannels==Channels(i))
        if find(MarkedChannels==Channels(i))
            h1=plot(1.8,-70,'bd-');
            set(h1,'MarkerSize',10);
            set(h1,'MarkerFaceColor','green');
        end
        if find(MarkedChannels2==Channels(i))
            a=find(MarkedChannels2==Channels(i));
            %h1=plot(1.3,-70,'bo-');
            %set(h1,'MarkerSize',10);
            %set(h1,'MarkerFaceColor','red');
            numerki=[1 3 4 5 6 2];
            numerki=[1:67];
            h=text(2.1,-90,num2str(numerki(a)));
            set(h,'FontSize',FontSize);
        end
        %disp('ble')
        %end
    end
%hold on;
end

y=Data;