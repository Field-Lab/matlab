clear
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes

MovieFilePath='I:\analysis\slices\2013-12-12-3-PH\movie005';
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
vreak
full_path='I:\analysis\slices\2013-12-12-3-PH\data005'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

g=importdata('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\NeuronsElectrodesv2.txt');

L=600;

for Event=65:length(g)
    NeuronID=g(Event,1);
    Electrode=g(Event,2);
    if find(idList==NeuronID) & Electrode<5
    %tutaj znalezc reording electrode
    SeedEl = neuronFile.getNeuronIDElectrode(NeuronID);
    Electrodes=electrodeMap.getAdjacentsTo(SeedEl,1)';
    
    spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times
    NeuronEI=NS512_LoadEIFromRawData(full_path,spikeTimes(10:min(509,length(spikeTimes))),60,20);
    
    EIforElectrodes=NeuronEI(Electrodes,:);
    PP=max(EIforElectrodes')-min(EIforElectrodes')
    RecordingElectrodeIndex=find(PP==max(PP))
    RecordingElectrode=Electrodes(RecordingElectrodeIndex)
    
    MaxY=ceil(max(NeuronEI(RecordingElectrode,:))*2/10)*10
    MinY=floor(min(NeuronEI(RecordingElectrode,:))*1.5/10)*10
    
    figure(10)
    clf;
    
    fid=fopen(['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\dane_combined\ID=' num2str(NeuronID) 'c'],'r','ieee-le')
    a=fread(fid,'int32');
    l=length(a);
    b=reshape(a,5,l/5); %MovieNumber,RepetitionNumber,PatternNumber,Latency,TimeRelativeToRepetitionBegin
    fclose(fid);
    
    MoviesForThisNeuron=[[1:8]+(Electrode-1)*8+32]
    for movie=1:length(MoviesForThisNeuron)
        SpikesForThisMovie=find(b(1,:)==MoviesForThisNeuron(movie))
        MovieDataFull=NS_MovieData_GlobalPath(MovieFilePath,MoviesForThisNeuron(movie),NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieDataFull);
        SE=MovieData(2:3:length(MovieData));
        TimeOfLastPulse=MovieData(length(MovieData)-2);
        
        spikes=zeros(50,L);
        for repetition=1:50
            PulseTime=MovieBegin+(repetition-1)*RepetPeriod+TimeOfLastPulse;
            d0=rawFile.getData(PulseTime,L)';
            d1=d0(RecordingElectrode+1,:);
            spikes(repetition,:)=d1;
        end
        subplot(4,2,movie);
        h=plot(spikes');
        set(h,'Color','b');        
        hold on
        grid on
        h11=gca;
        set(h11,'YLim',[MinY MaxY]);
        
        hj=0;
        for repetition=1:50
            SpikesForThisRepetition=find(b(2,:)==repetition);
            SpikesForThisMovieAndRepetition=intersect(SpikesForThisMovie,SpikesForThisRepetition)
            hj=hj+length(SpikesForThisMovieAndRepetition)
            for SuspectedSpike=1:length(SpikesForThisMovieAndRepetition)
            
                TimeOfSpike=b(5,SpikesForThisMovieAndRepetition(SuspectedSpike));
                LatencyOfSpike=TimeOfSpike-TimeOfLastPulse
            
                if LatencyOfSpike>0 & LatencyOfSpike<L-200
                    PulseTime=MovieBegin+(repetition-1)*RepetPeriod+TimeOfLastPulse
                    d0=rawFile.getData(PulseTime,L)';
                    d1=d0(RecordingElectrode+1,:);
                    h1=plot([LatencyOfSpike-10:LatencyOfSpike+20],d1([LatencyOfSpike-10:LatencyOfSpike+20]));
                    %h1=plot(d1)
                    set(h1,'Color','r');
                end
            end
        end                                                                
    end
    FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\report_2014_02_13\figures2\';
    FullImageName=[FigurePath 'Neuron' num2str(NeuronID) '_Stim' num2str(Electrode) 'add.tif'];
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 9]);
    set(h,'PaperPosition',[0 0 16 9]); 
    print(h, '-dtiff', '-r120', FullImageName);
    end
end    

%NeuronEI=NS512_LoadEIFromRawData(RawDataPath,SpikeTimesForElectrode,Offset
%+60,Offset);