clear
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
MovieFile='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\MovieFiles\2013-12-12-3\movie001';
[NumberOfMovies,NumberOfPatternsPerMovie,AllPatterns,Patterns]=NS512_MoviePatterns(MovieFile,NS_GlobalConstants);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
SPFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\sp_files-2015-05-24\';

%clear qw;
PatternsToAnalyze=unique(AllPatterns)';

AllSpikes=zeros(25,512);
AllSpikes2=zeros(25,length(PatternsToAnalyze),512);
%ManySpikes=zeros(25,512);
for Amplitude=1:25
    Amplitude
    SpikesForAmplitude=zeros(1,512);
    ManySpikesForAmplitude=zeros(1,512);
    for i=1:length(PatternsToAnalyze)
        Pattern=PatternsToAnalyze(i);
        MoviesForPattern=find(Patterns(1:450,Pattern)>0);
        ChannelsToExclude=electrodeMap.getAdjacentsTo(Pattern,3)';    
    
    %for Amplitude=2:25        
        Movie=MoviesForPattern(Amplitude);
        fid=fopen([SPFilesPath 'sp_p' num2str(Pattern) 'm' num2str(Movie)],'r','ieee-le');
        a=fread(fid,'int32');
        b=reshape(a,length(a)/3,3);
        fclose(fid);  
        
        Spikes=zeros(1,512);
        for el=1:512
            l=find(b(:,1)==el);
            if l
                SpikeTimesHistogram=hist(b(l,3),[0.5:1:599.5]);
                TimeLocking=NS512_AreTheSpikesTimeLocked_NoSpikeNumberLimit(SpikeTimesHistogram,40);
                if TimeLocking==1
                    Spikes(el)=length(l);     
                else
                    Spikes(el)=-length(l); 
                end
                
                %AllSpikes2(Amplitude,i,el)=length(l);
            end
        end
        Spikes(ChannelsToExclude)=0;
        AllSpikes2(Amplitude,i,:)=Spikes;
        Spikes2=find(Spikes>25);
        SpikesForAmplitude=SpikesForAmplitude+Spikes;
        %SpikesOverall=max(SpikesOverall,Spikes);
        %Amplitude        
    end    
    AllSpikes(Amplitude,:)=AllSpikes(Amplitude,:)+SpikesForAmplitude;
end