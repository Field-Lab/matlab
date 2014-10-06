function [AmplitudesVsChannels,RMSVsChannelsAll,SamplesOutput]=NS512_FindThresholdForLocalStimulation4AK_PH2012_04_29(DataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Samples,SamplesNumber);
% Only one pattern number, but range of movies and also specified group of
% electrodes to look for spikes on
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

ChannelsToRead=Channels;
AmplitudesVsChannels=[];
RMSVsChannelsAll = [];

NumberOfSpikesInFile=2400000;
SamplesOutput=SamplesNumber;
ArtifactThresholdNumber=5;
SpikesNumberThreshold=25;

N=7;
for MovieNumber=Movies
    mn=MovieNumber;    
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);   %read data                 
    sdfsdf=size(DataTraces0);
    DataTraces=DataTraces0(1:45,Channels,Samples); %take only some samples from the data
    sdfsdf2=size(DataTraces);
      
    % 1. Artifact estimation
    [Artifact]=NS512_ArtifactEstimationSomaticManyChannels(DataTraces,ArtifactThresholdNumber); %size of Artifact is 1 x 512 x 70
    % 2. Artifact subtraction 
    TracesWithoutArtifact=NS512_SubtractArtifact(DataTraces,Artifact); %size of TracesWithoutArtifact is 100 x 512 x 7
    % 3. Spike detection  
    
    [Events,LatencyCriterion,HalfAmpCrossingIndexes]=NS512_FindSpikes(TracesWithoutArtifact,ChannelsToRead,-18,3,[3 20],[2 40]);
    % 4. Identify channels with high response rate
    ChannelsWithSpikes=NS512_FindChannelsWithHighResponseRate(Events,SpikesNumberThreshold);
    RealChannelsWithSpikes=ChannelsToRead(ChannelsWithSpikes);
    SpikesTimings=HalfAmpCrossingIndexes(:,ChannelsWithSpikes);      
                    
    for i=1:length(ChannelsWithSpikes)
        PatternNumber
        mn
        ChannelsWithSpikes
        Channel=ChannelsToRead(ChannelsWithSpikes(i)); %ChannelsToRead - to samo co wejsciowe Channels, typowo wszystkie kana?y matrycy
        ChannelsPlot = electrodeMap.getAdjacentsTo(Channel,1); % Elektrody sasiadujace z elektroda Channel (ta z wiecej niz 50 spikami)        
        
        WaveformTypes=Events(:,ChannelsWithSpikes(i)); %konkretny kana?, jeden z tych z du?? ilo?ci? spików
        TracesWithSpikesInd=find(WaveformTypes==1);
        size(TracesWithSpikesInd);
        ST=SpikesTimings(:,i);        
        [CorrectedTraces,EI,UniSpikesIndicCorrected]=NS512_TimingsForDetectedNeuron2(TracesWithoutArtifact(TracesWithSpikesInd,:,:),ChannelsPlot); % TracesWithoutArtifact - all the traces, not only with spikes!                      
        quality=NS512_PlotDetectedSpikes(TracesWithoutArtifact,ChannelsPlot,WaveformTypes,CorrectedTraces,EI);  %zmiana     
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
        
        %To co dalej nie zale?y od tego, czy mamy dobr?jako??!!! pewnie
        %docelowo zmieni?
                
         TraceSS=TracesWithoutArtifact(TracesWithSpikesInd,:,:); % wiersze z zarejestrowanymi spikami
  
         samples512x35xN=[];
         parametry=[];
        for j=1:size(TracesWithSpikesInd)
            samples512x35=[];
            for k=1:length(Samples)
                samples512x1=reshape(TraceSS(j,:,k),512,1);                                               
                
                SamplesOutput=SamplesOutput+1; %ilosc probek++
                samples512x35=cat(2,samples512x35,samples512x1);
                %jeœli przekroczymy max ilosc probek
                if mod(SamplesOutput,NumberOfSpikesInFile)==0        %maksymalna ilosc probek to 2 400 000                    
                    nazwa = ['data00000'  num2str(floor((SamplesOutput-1)/NumberOfSpikesInFile))  '.bin']
                    fid=fopen(nazwa,'a');
                    samples512x35xN=cat(2,samples512x35xN,samples512x35);
                    samples512x35xN16bit=int16(samples512x35xN); %rzutowanie
                    size(samples512x35xN16bit);
                    samples512x35xN12bit = NS512_CompressData16bit12bit(cat(1, zeros(1,size(samples512x35xN16bit,2)),samples512x35xN16bit)); %kompresja 513x35xN
                              
                    fwrite(fid, samples512x35xN12bit,'uint8'); %przerobic na uint8!!!!!!
                    clear samples512x35xN16bit;
                    fclose(fid); 
                
                    samples512x35=[];
                    samples512x35xN=[];                  
                end                
            end      
            parametry=cat(1,parametry,[PatternNumber,mn,TracesWithSpikesInd(j),Channel]); %dodac numer elekrtody ze spikiem
            samples512x35xN=cat(2,samples512x35xN,samples512x35);               
        end                  
            nazwa = ['data00000'  num2str(floor(SamplesOutput/NumberOfSpikesInFile))  '.bin'];  
              
            fid=fopen(nazwa,'a');
            samples512x35xN=int16(samples512x35xN);
            samples512x35xN12bit = NS512_CompressData16bit12bit(cat(1, zeros(1,size(samples512x35xN,2)),samples512x35xN)); %kompresja 513x35xN z dopisaniem wiersza zer
            fwrite(fid, samples512x35xN12bit,'uint8'); %przerobic na uint8!!!!!!
            fclose(fid);
                 
            fid=fopen('parametry.bin','a');
            fwrite(fid,parametry','uint16'); %odczyt: f=fopen('parametry.bin','r'); c=fread(f,'uint16'); N=length(c); d=reshape(c,4,N/4)';
            fclose(fid);                
    end
    ChannelsToRead=NS_RemoveBadChannels(ChannelsToRead,RealChannelsWithSpikes);
end