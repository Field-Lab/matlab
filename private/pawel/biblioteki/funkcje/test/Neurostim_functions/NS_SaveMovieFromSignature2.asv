function M=NS_SaveMovieFromSignature2(Traces,Channels,ElectrodeToShow,OversamplinFactor,ArrayID,WritePathFigs,NS_GlobalConstants);
%Traces - array NxT, where N - number of channels, T - number of samples in
%each traces
%ElectrodeToShow=37;

STraces=size(Traces);
ChannelsNumber=length(Channels);
NumberOfFrames=STraces(2);

t=1:NumberOfFrames;
t2=1:1/OversamplinFactor:NumberOfFrames;
TracesNew=interp1(t,Traces',t2,'spline')';
clear Traces;
Traces=TracesNew;
STraces=size(Traces);
NumberOfFrames=STraces(2);

for i=1:ChannelsNumber
    Traces(i,:)=Traces(i,:)-mean(Traces(i,1:4));
end

FigureNumber=100;

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

figure(FigureNumber);
clf;
maxes=max(max(abs(Traces')));
for frame=1:NumberOfFrames*OversamplinFactor %for each frame...
    clf;
    for i=1:ChannelsNumber
        Kanal=Channels(i)
        
        X=(Coordinates(1,i)-X0)/Xspread+0.05;
        Y=(Coordinates(2,i)-Y0)/Yspread+0.1;
        signal=Traces(i,frame)/maxes*150;%/maxes(i)*20;    
        %h=plot(X,Y,'bo');
        %text(X,Y,num2str(Channels(i)));
        hold on;
        if Channels(i)==ElectrodeToShow
            h=text(X,Y,num2str(ElectrodeToShow));
            set(h,'FontSize',16);
        end
        if i==7
            %text(X,Y,'7');
        end
        rad=round(abs(signal)/3);
        %rad=round(sqrt(abs(signal)/4))^2;

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
    h=text(0,1,['t=' num2str(frame/(20*OversamplinFactor),'%6.2f')]);
    set(h,'FontSize',16);
    axis([-0.1 1.1 -0.1 1.1]);
    h=gca;
    set(h,'Position',[0.05 0.05 0.9 0.9]);
    set(h,'Visible','Off');
        o
    h=gcf;
    FullName=[WritePathFigs '\' 'frame' num2str(frame)];            
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[9 10]);
    set(h,'PaperPosition',[0 0 9 10]); 
    print(h, '-dtiff', '-r80', FullName);
    
    %M(frame)=getframe;
    refresh;
    pause(0.1);       
end
movie2avi(M,'proba','FPS',10);