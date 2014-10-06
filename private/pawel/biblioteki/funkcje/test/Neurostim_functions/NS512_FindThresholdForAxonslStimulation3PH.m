function [AmplitudesVsChannels,RMSVsChannelsAll]=NS512_FindThresholdForAxonalStimulation3PH(DataPath,ArtifactDataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Samples);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
SpikesNumberThreshold=40;

ChannelsToRead=Channels;
AmplitudesVsChannels=[];
RMSVsChannelsAll = [];

ArtifactThresholdNumber=10;
SpikesNumberThreshold=30;

N=7;
for MovieNumber=Movies
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,1,PatternNumber,MovieNumber,0,0);                
    DataTraces=DataTraces0(1:100,[1:512],Samples);
    TracesWithoutArtifact=DataTraces; %because we use TTx based artifact recording!
    
    % 1. Artifact estimation
    %[Artifact]=NS512_ArtifactEstimationSomatic(DataTraces,ArtifactThresholdNumber); %size of Artifact is 1 x 512 x 70
    % 2. Artifact subtraction 
    %TracesWithoutArtifact=NS512_SubtractArtifact(DataTraces,Artifact); %size of TracesWithoutArtifact is 100 x 512 x 7
    % 3. Spike detection
    [Events,LatencyCriterion,HalfAmpCrossingIndexes]=NS512_FindSpikes(TracesWithoutArtifact,ChannelsToRead,-25,3,[5 20],[2 40]);
    % 4. Identify channels with high response rate
    ChannelsWithSpikes=NS512_FindChannelsWithHighResponseRate(Events,SpikesNumberThreshold);
    RealChannelsWithSpikes=ChannelsToRead(ChannelsWithSpikes);
    SpikesTimings=HalfAmpCrossingIndexes(:,ChannelsWithSpikes);
                
    for i=1:0 %length(ChannelsWithSpikes)
        Channel=ChannelsToRead(ChannelsWithSpikes(i));                
        ChannelsPlot = electrodeMap.getAdjacentsTo(Channel,1); % Elektrody sasiadujace z elektroda Channel (ta z wiecej niz 50 spikami)        
        
        WaveformTypes=Events(:,ChannelsWithSpikes(i));
        TracesWithSpikesInd=find(WaveformTypes==1);
        ST=SpikesTimings(:,i);
        [CorrectedTraces,EI,UniSpikesIndicCorrected]=NS512_TimingsForDetectedNeuron2(TracesWithoutArtifact(TracesWithSpikesInd,:,:),ChannelsPlot); % TracesWithoutArtifact - all the traces, not only with spikes!
        quality=NS512_PlotDetectedSpikes(TracesWithoutArtifact,ChannelsPlot,WaveformTypes,CorrectedTraces,EI);       
        if quality==1
            WritePathFigs=WritePathFigsGood;
        else
            WritePathFigs=WritePathFigsBad;
        end
        h=gcf;
        FullName=[WritePathFigs '\' 'p' num2str(PatternNumber) '_m' num2str(MovieNumber) '_el' num2str(Channel)];            
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]); 
        print(h, '-dtiff', '-r120', FullName);
    end
    ChannelsToRead=NS_RemoveBadChannels(ChannelsToRead,RealChannelsWithSpikes);
end