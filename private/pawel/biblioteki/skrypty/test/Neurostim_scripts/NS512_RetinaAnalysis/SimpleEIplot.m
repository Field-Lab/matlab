DataPath='G:\D\Home\Pawel\analysis\retina\2012sept\2012-09-24-0\scan_proc'

PatternNumber=46;
MovieNumber=63;

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);
Samples=[16:55];
EI=reshape(mean(DataTraces(1:50,:,16:55)),512,length(Samples));
figure(1)
clf
for i=1:512
    EI(i,:)=EI(i,:)-mean(EI(i,length(Samples)-4:length(Samples)));
end
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
ElectrodesToExclude=electrodeMap.getAdjacentsTo(PatternNumber,1)'

EI(ElectrodesToExclude,:)=5;
EI(PatternNumber,:)=150
h=NS512_ShowEIAsCircles(EI/2,500,[1:512],PatternNumber,[-1000 1000],[-500 500]);