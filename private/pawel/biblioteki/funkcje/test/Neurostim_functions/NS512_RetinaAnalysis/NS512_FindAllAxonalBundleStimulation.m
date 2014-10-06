function Thresholds=NS512_FindAllAxonalBundleStimulation(DataPath,Patterns,Movies,Threshold,SamplesToAnalyze)
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
%Threshold - wartosc ujemna, bo szukamy negatywnej wartosci
Thresholds=zeros(length(Patterns),512);

for i=1:length(Patterns)
    ElectrodesToExclude=electrodeMap.getAdjacentsTo(Patterns(i),1)';
    ThresholdsForPattern=zeros(1,512);
    for j=1:length(Movies)
        [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Patterns(i),Movies(j),0,0);
        DataTraces=DataTraces0(1:50,[1:512],SamplesToAnalyze);
        d=reshape(mean(DataTraces),512,length(SamplesToAnalyze));   
        means=mean(d(:,length(SamplesToAnalyze)-10:length(SamplesToAnalyze)),2);
        for k=1:512
            d(k,:)=d(k,:)-means(k);
        end        
        PP=min(d');
        
        ElectrodesAboveThreshold=PP<-Threshold;
        ElectrodesAboveThreshold(ElectrodesToExclude)=0;
        ElectrodesWithoutThreshold=ThresholdsForPattern==0;
        ElectrodesToUpdate=find(ElectrodesAboveThreshold.*ElectrodesWithoutThreshold==1);                                
        ThresholdsForPattern(ElectrodesToUpdate)=j;        
    end
    Thresholds(i,:)=ThresholdsForPattern;
end        