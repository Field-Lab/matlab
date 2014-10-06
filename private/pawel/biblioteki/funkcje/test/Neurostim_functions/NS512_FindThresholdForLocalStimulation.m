function MovieNumber=NS512_FindThresholdForLocalStimulation(DataPath,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels);

for MovieNumber=Movies    
    MovieNumber
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);    
    size(DataTraces0)
    DataTraces=DataTraces0(1:100,Channels,1:50);
    [Events,Artifact]=NS512FindResponsesToLocalStimulation(DataTraces,25,3,10); %which traces include spikes that exceed threshold (independently for each channel)
    EventsPerChannel=sum(Events); %how many spikes on each electrode
    a=max(EventsPerChannel);
    if a>50 %if there is electrode with more than 50 spikes...
        [p,index]=find(EventsPerChannel==a) %which electrodeis it?
        WaveformTypes=Events(:,index(1)); %0 - artifact only 1 - including spike
        artifacts=find(WaveformTypes==0); %find all the "artifact only" traces
        Artifact=mean(DataTraces(artifacts,:,:)); % find EI of the artifact
        Traces=NS512_SubtractArtifact(DataTraces,Artifact); 
        y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(Traces,Channels,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels);        
        break;
    end
end
if a<51
    MovieNumber=0;
end