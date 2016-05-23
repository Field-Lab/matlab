NS_GlobalConstants=NS_GenerateGlobalConstants(500);
DataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
Movies=[1:2:63];
FigurePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\summary4';

%ThresholdMovies=zeros(1,PatternsToShow);
Threshold=80;

figure(3)
clf

AmplitudesLabels={};
for i=1:length(Movies)    
    StimAmplitudes(i)=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],Patterns(1),Movies(i),NS_GlobalConstants);
    if i/2==round(i/2)
        AmplitudesLabels{i}=num2str(StimAmplitudes(i),'%4.2f');
    else
        AmplitudesLabels{i}='';
    end
end
break
AmplitudeMin=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],Patterns(1),Movies(1),NS_GlobalConstants)
AmplitudeMax=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],Patterns(1),max(Movies),NS_GlobalConstants)
SamplesToAnalyze=[7:136];
t=[1:130]/20;
for i=1:length(Patterns)
    Pattern=Patterns(i);
    ElectrodesToExclude=electrodeMap.getAdjacentsTo(Pattern,1)';
    %step 1: which electrodes to consider - 
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Patterns(i),max(Movies),0,0);
    DataTraces=DataTraces0(1:50,[1:512],SamplesToAnalyze);
    d=reshape(mean(DataTraces),512,130);   
    means=mean(d1(:,111:130),2);
    for k=1:512
        d(k,:)=d(k,:)-means(k);
    end
    PP=min(d');     
    ElectrodesAboveThreshold=PP<-Threshold
    ElectrodesAboveThreshold(ElectrodesToExclude)=0
    ElectrodesToAnalyze=find(ElectrodesAboveThreshold==1)
    for el=1:length(ElectrodesToAnalyze)
        clf
        subplot(2,2,1);
        h=NS512_ShowEIFrameAsCircles(PP'/15,500,[1:512],Patterns(i),ElectrodesToAnalyze(el),[],[-1005 1005],[-505 505]);
        subplot(2,2,2);
        plot(t,d(ElectrodesToAnalyze(el),:));        
   
        for j=1:length(Movies)
            Amplitude=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],Patterns(i),Movies(j),NS_GlobalConstants);
            [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Patterns(i),Movies(j),0,0);
            DataTraces=DataTraces0(1:50,[1:512],SamplesToAnalyze);
            d=reshape(mean(DataTraces),512,130);    
            for k=1:512
                d(k,:)=d(k,:)-means(k);
            end
            PP=min(d'); 
            %PP=max(d')-min(d'); 
            A1(j)=abs(PP(ElectrodesToAnalyze(el)));
        end
        subplot(2,2,3);
        loglog(StimAmplitudes,A1,'bd-');
        h1=gca;
        set(h1,'XLim',[AmplitudeMin AmplitudeMax]);
        set(h1,'YLim',[5 2000]);
        set(h1,'XTick',StimAmplitudes);
        set(h1,'XTickLabels',AmplitudesLabels);
        grid on;
        FigureName=['C:\pawel\nauka\analiza\retina\2012-09-27-4\summary4\stim' num2str(Pattern) '_rec' num2str(ElectrodesToAnalyze(el))]                
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]); 
        print(h, '-dtiff', '-r100', FigureName);  
    end
end