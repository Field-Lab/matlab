function y=NS_PlotManySignaturesOnArrayLayout_SP(Traces,Channels,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants);
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

Fs=NS_GlobalConstants.SamplingFrequency;

FigureNumber=FigureProperties.FigureNumber;
TimeRange=FigureProperties.TimeRange;
AmplitudeRange=FigureProperties.AmplitudeRange;
FontSize=FigureProperties.FontSize;
Colors=FigureProperties.Colors;
LineWidth=FigureProperties.LineWidth;

FontSize=20;

%XTicks=[0 0.5 1 1.5 2 2.5 3];
%YTicks=[-200 -100 0 100 200];
%if AmplitudeRange(2)==100
%    YTicks=[-150 -100 -50 0 50 100];
%end

Xstep=60; 
Ystep=60;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
TimeRange=TimeRange./Fs*1000;

t=[TimeRange(1):(1/Fs*1000):TimeRange(2)];
%SSignature=size(Signature);
t=([0:STraces(3)-1])*1000/Fs;
%TimeRange(1)
%(1-TimeRange(1))*1000/Fs
%t(1)

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

%figure(FigureNumber);
%clf;

for k=1:STraces(1) %for each signature...
    k;
    for i=1:ChannelsNumber
        %iteration=k
        i;
        X=(Coordinates(1,i)-X0)/Xspread+0.15; %0.05;
        Y=(Coordinates(2,i)-Y0)/Yspread+0.12; %0.1;
        XSize=Xstep/Xspread;
        YSize=Ystep/Yspread;
        %subplot('Position',[X+(i-1)*0.03-0.01 Y XSize*0.9 YSize*0.75]);
        signal0=Traces(k,i,:);       
        signal=reshape(signal0,1,STraces(3));  
        WaveformColorIndex=WaveformTypes(k)-floor(WaveformTypes(k)/length(Colors))*length(Colors)+1;   
        %kolor=Colors(WaveformColorIndex);
        if (Colors(WaveformColorIndex)~='n')         
            a=plot(t,signal,Colors(WaveformColorIndex));
            if Colors(WaveformColorIndex)=='g'
                set(a,'Color',[0 0.75 0]); 
            end
            
            set(a,'LineWidth',LineWidth);
            %pause(1)
            %refresh(FigureNumber)
            %set(a,'LineWidth',2);
            hold on;        
            %if k==STraces(1)
            %  h=gca;
            %  AG=get(h,'YLim');
            %  text(TimeRange(1)+(TimeRange(2)-TimeRange(1))*0.85,AG(1)+(AG(2)-AG(1))*0.9,num2str(Channels(i)));
            %  %text(TimeRange(1)+(TimeRange(2)-TimeRange(1))*0.8,AmplitudeRange(1)+(AmplitudeRange(2)-AmplitudeRange(1))*0.8,num2str(Channels(i)));
            %end
            xlim(TimeRange);
            ylim(AmplitudeRange);
            grid on;
            xt=TimeRange(1)+(TimeRange(2)-TimeRange(1))*0.8;
            yt=AmplitudeRange(1)+(AmplitudeRange(2)-AmplitudeRange(1))*0.9;
            %h=text(2.2,150,num2str(Channels(i)));
            %set(h,'FontSize',FontSize);
            h=gca;  
            set(h,'XTick',[0:0.5:2.5]);
            set(h,'YTick',[-250:50:250]);
            set(h,'YTickLabel',{'' '-200' '' '-100' '' '0' '' '100' '' '200' ''});
            
            %set(h,'YTick',[-150:25:100]);
            %set(h,'YTickLabel',{'-150' '' '-100' '' '-50' '' '0' '' '150' '' '100' ''});
            
            %set(h,'YTick',YTicks);
            if i==SE        
                set(h,'FontSize',FontSize);                
                h=xlabel('Time [ms]');
                set(h,'FontSize',FontSize);
                h=ylabel('Amplitude [\muV]');
                set(h,'FontSize',FontSize);
            else
                set(h,'FontSize',FontSize);
                %set(h,'XTickLabel',[0.5 1 1.5 2 2.5]);
                %set(h,'XTickLabel',[]);
                set(h,'YTickLabel',[]);
                h=xlabel('Time [ms]');
                set(h,'FontSize',FontSize);
            end
        end
        if k==STraces(1)
              h=gca;
              AG=get(h,'YLim');
              %text(TimeRange(1)+(TimeRange(2)-TimeRange(1))*0.85,AG(1)+(AG(2)-AG(1))*0.9,num2str(Channels(i)));
              %text(TimeRange(1)+(TimeRange(2)-TimeRange(1))*0.8,AmplitudeRange(1)+(AmplitudeRange(2)-AmplitudeRange(1))*0.8,num2str(Channels(i)));
        end
    end
    %refresh;
    %pause(1);
end

y=Coordinates;