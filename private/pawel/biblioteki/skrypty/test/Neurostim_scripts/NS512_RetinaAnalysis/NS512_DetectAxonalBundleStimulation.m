NS_GlobalConstants=NS_GenerateGlobalConstants(500);
DataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
Movies=[1:2:63];

ThresholdMovies=zeros(1,PatternsToShow);
Threshold=20;

figure(3)
clf

for i=1:length(Movies)
    StimAmplitudes(i)=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],Patterns(1),Movies(i),NS_GlobalConstants);
end

SamplesToAnalyze=[8:137];
%ChoosenElectrode=473;
for i=1:1%length(Patterns)
    Pattern=Patterns(i);
    ElectrodesToExclude=electrodeMap.getAdjacentsTo(Pattern,1)';
    %step 1: which electrodes to consider - 
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Patterns(i),max(Movies),0,0);
    DataTraces=DataTraces0(1:50,[1:512],SamplesToAnalyze);
    d=reshape(mean(DataTraces),512,130);    
    PP=max(d')-min(d'); 
    
    ElectrodesAboveThreshold=PP>Threshold;
    ElectrodesAboveThreshold(ElectrodesToExclude)=0;
    ElectrodesToAnalyze=find(ElectrodesMarkedToPlot==1)
    for el=1:length(ElectrodesToAnalyze)
        subplot(2,1,1);
        h=NS512_ShowEIFrameAsCircles(PP'/5,500,[1:512],Patterns(i),ElectrodesToAnalyze(i),[],[-1005 1005],[-505 505]);
   
        for j=1:length(Movies)
            Amplitude=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],Patterns(i),Movies(j),NS_GlobalConstants)
            [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Patterns(i),Movies(j),0,0);
            DataTraces=DataTraces0(1:50,[1:512],SamplesToAnalyze);
            d=reshape(mean(DataTraces),512,130);    
            PP=max(d')-min(d'); 
            A1(j)=PP(402);
            ElectrodesMarkedToPlot=PP>Threshold;
            ElectrodesMarkedToPlot(ElectrodesToExclude)=0;
            ElectrodesToPlot=find(ElectrodesMarkedToPlot==1)
            subplot(2,1,2);
            if ElectrodesToPlot
                h=plot(Amplitude,PP(ElectrodesToPlot),'bd');
                hold on
            end 
        end
    end
end
AmplitudeMin=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],Patterns(1),Movies(1),NS_GlobalConstants)
AmplitudeMax=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],Patterns(1),max(Movies),NS_GlobalConstants)
grid on
h1=gca
set(h1,'XLim',[AmplitudeMin AmplitudeMax]);
set(h1,'XScale','log')
set(h1,'YScale','log')
set(h1,'XTick',StimAmplitudes)