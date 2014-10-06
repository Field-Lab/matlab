function [AmplitudesVsChannels,RMSVsChannelsAll]=NS512_FindThresholdForLocalStimulation4PH(DataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Samples);
% Only one pattern number, but range of movies and also specified group of
% electrodes to look for spikes on
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

ChannelsToRead=Channels;
AmplitudesVsChannels=[];
RMSVsChannelsAll = [];

ArtifactThresholdNumber=5;
SpikesNumberThreshold=25;

N=7;
for MovieNumber=Movies
    mn=MovieNumber
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);   %read data                 
    sdfsdf=size(DataTraces0)
    figure(111)
    plot(ChannelsToRead,'bd-')
    figure(112)
    plot(Channels,'bd-')
    DataTraces=DataTraces0(1:45,Channels,Samples); %take only some samples from the data
    sdfsdf2=size(DataTraces)
    
    % 1. Artifact estimation
    [Artifact]=NS512_ArtifactEstimationSomatic(DataTraces,ArtifactThresholdNumber); %size of Artifact is 1 x 512 x 70
    % 2. Artifact subtraction 
    TracesWithoutArtifact=NS512_SubtractArtifact(DataTraces,Artifact); %size of TracesWithoutArtifact is 100 x 512 x 7
    figure(11)
    % 3. Spike detection
    size(DataTraces)
    ChannelsToRead
    size(TracesWithoutArtifact)
    [Events,LatencyCriterion,HalfAmpCrossingIndexes]=NS512_FindSpikes(TracesWithoutArtifact,ChannelsToRead,-18,3,[3 20],[2 40]);
    % 4. Identify channels with high response rate
    ChannelsWithSpikes=NS512_FindChannelsWithHighResponseRate(Events,SpikesNumberThreshold);
    RealChannelsWithSpikes=ChannelsToRead(ChannelsWithSpikes);
    SpikesTimings=HalfAmpCrossingIndexes(:,ChannelsWithSpikes);
                
    for i=1:length(ChannelsWithSpikes)
        Channel=ChannelsToRead(ChannelsWithSpikes(i));                
        ChannelsPlot = electrodeMap.getAdjacentsTo(Channel,1); % Elektrody sasiadujace z elektroda Channel (ta z wiecej niz 50 spikami)        
        
        WaveformTypes=Events(:,ChannelsWithSpikes(i))
        size(WaveformTypes)
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
        set(h,'PaperSize',[13 9]);
        set(h,'PaperPosition',[0 0 13 9]); 
        print(h, '-dtiff', '-r120', FullName);
    end
    ChannelsToRead=NS_RemoveBadChannels(ChannelsToRead,RealChannelsWithSpikes);
end