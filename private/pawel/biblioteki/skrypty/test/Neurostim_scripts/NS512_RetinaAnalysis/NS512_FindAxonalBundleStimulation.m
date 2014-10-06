clear
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
ArrayID=500;
DataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
Movies=[1:2:63];

ThresholdFilePath1=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\stim_scan\thresholds_1']
fid1=fopen(ThresholdFilePath1,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons1=NeuronInformation(:,1);
Electrodes1=NeuronInformation(:,3);
Amplitudes1=NeuronInformation(:,4);

ThresholdFilePath2=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\stim_scan\thresholds_2']
fid1=fopen(ThresholdFilePath2,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons2=NeuronInformation(:,1);
Electrodes2=NeuronInformation(:,3);
Amplitudes2=NeuronInformation(:,4);

Neurons=[Neurons1' Neurons2']';
Electrodes=[Electrodes1' Electrodes2']';
Amplitudes=[Amplitudes1' Amplitudes2']';

Patterns=sort(Electrodes);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
Channels=[1:512];
Coordinates=zeros(2,length(Channels));
ChannelsNumber=length(Channels);
for i=1:ChannelsNumber
    Coordinates(1,i)=electrodeMap.getXPosition(Channels(i));
    Coordinates(2,i)=electrodeMap.getYPosition(Channels(i));
end

N=length(Patterns)
break
%OutputFileName='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\proby2\proba.gif';
for i=[25 58 64]%[1:6 8:11 13:22 24:26 28:32 34:54 56:70]
    figure(4)
    clf
    OutputFileName1=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\proby5\p' num2str(i) '.gif'];
    %ElectrodesToLookAt=find(Coordinates(2,:)<Coordinates(2,Patterns(i))-120);
    for f=1:0%31
        % * * * * Dodac kodowanie szerokosci sygna³u kolorem oraz informacjê o
        % amplitudzie pr¹du
        clf
        [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Patterns(i),(f-1)*2+1,0,0);
        Amplitude=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],Patterns(i),(f-1)*2+1,NS_GlobalConstants)
        DataTraces=DataTraces0(1:50,[1:512],[8:137]);
        d=reshape(mean(DataTraces),512,130);
    
        PP=max(d')-min(d'); 
        %length(find(max(PP(ElectrodesToLookAt)>50)))
        h=NS512_ShowEIFrameAsCircles(PP'/5,500,[1:512],Patterns(i),[],[],[-1005 1005],[-505 505]);
        set(h,'Visible','off');
        h1=text(-900,450,['I=' num2str(Amplitude,'%8.2f')]);
        set(h1,'FontSize',20);
        frame = getframe(h);
        im = frame2im(frame);
        [imind,map] = rgb2ind(im,256);
        if f == 1
            imwrite(imind,map,OutputFileName1,'gif', 'Loopcount',1);
        else
            imwrite(imind,map,OutputFileName1,'gif','WriteMode','append','DelayTime',0.2);
        end
        %pause(1)
        %refresh
        %max(PP(ElectrodesToLookAt));
    end
    
    figure(5)
    clf
    %OutputFileName2=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\proby5\pa' num2str(i) '.gif'];
    OutputFileName2=['C:\pawel\nauka\analiza\retina\2012-09-27-4\summary3\figures_specific\pa' num2str(i) '.gif'];
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Patterns(i),63,0,0);
    DataTraces=DataTraces0(1:50,[1:512],[1:140]);
    d=reshape(mean(DataTraces),512,140);
    SEI=size(d);
    for j=1:SEI(1)
            d(j,:)= d(j,:)-mean([d(j,1) d(j,SEI(2))]);
        end
    PP=max(d')-min(d');
    for f=1:45%135
        clf        
        h=NS512_ShowEIFrameAsCircles(d(:,f+5)/2,500,[1:512],Patterns(i),[],[],[-1005 1005],[-505 505]);
        set(h,'Visible','off');
        h1=text(-900,450,['t=' num2str((f+5)/20,'%8.2f')]);
        set(h1,'FontSize',20);
        frame = getframe(h);
        im = frame2im(frame);
        [imind,map] = rgb2ind(im,256);
        if f == 1
            imwrite(imind,map,OutputFileName2,'gif', 'Loopcount',1);
        else
            imwrite(imind,map,OutputFileName2,'gif','WriteMode','append','DelayTime',0.1);
        end
    end
end