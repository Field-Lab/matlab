function M=NS_SaveMovieFromSignature(Traces,Channels,ElectrodeToShow,ArrayID,FigureProperties,NS_GlobalConstants);
%Traces - array NxT, where N - number of channels, T - number of samples in
%each traces
%ElectrodeToShow=37;

STraces=size(Traces);
ChannelsNumber=length(Channels);
NumberOfFrames=STraces(2);

for i=1:ChannelsNumber
    Traces(i,:)=Traces(i,:)-mean(Traces(i,1:4));
end

Fs=NS_GlobalConstants.SamplingFrequency;
FigureNumber=FigureProperties.FigureNumber;
TimeRange=FigureProperties.TimeRange;
AmplitudeRange=FigureProperties.AmplitudeRange;
FontSize=FigureProperties.FontSize;
Colors=FigureProperties.Colors;

Xstep=60; 
Ystep=60;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
%TimeRange=TimeRange./Fs*1000;
%t=[TimeRange(1):(1/Fs*1000):TimeRange(2)];
%SSignature=size(Signature);
%t=([0:STraces(3)-1])*1000/Fs;

Coordinates=zeros(2,length(Channels));
%ChannelsNumber=length(Channels);
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

figure(111);
plot(Traces')

figure(FigureNumber);
clf;
maxes=max(max(abs(Traces')));
for frame=1:NumberOfFrames %for each frame...
    clf;
    for i=1:ChannelsNumber
        X=(Coordinates(1,i)-X0)/Xspread+0.05;
        Y=(Coordinates(2,i)-Y0)/Yspread+0.1;
        signal=Traces(i,frame)/maxes*150;%/maxes(i)*20;    
        %h=plot(X,Y,'bo');
        %text(X,Y,num2str(Channels(i)));
        hold on;
        if i==ElectrodeToShow
            text(X,Y,num2str(ElectrodeToShow));
        end
        if i==7
            %text(X,Y,'7');
        end
        rad=round(abs(signal)/4);
        if rad~=0
            h=plot(X,Y,'bo');
            if signal<0
                cl=[1 0 0];
            else
                cl=[0 0 1];
            end      
            set(h,'MarkerEdgeColor',cl);
            set(h,'MarkerFaceColor',cl);
            set(h,'MarkerSize',rad); %abs(signal)*i);
        end
    end
    %hold on;
    text(0,1,['t=' num2str(frame/20,'%6.2f')]);
    axis([-0.2 1.2 -0.2 1.2]);
    M(frame)=getframe;
    refresh;
    pause(0.1);       
end
%movie2avi(M,'proba','FPS',10);