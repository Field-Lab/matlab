clear

MovieFilePath='I:\analysis\slices\2013-12-12-3-PH\movie005';
NS_GlobalConstants=NS_GenerateGlobalConstants(500);

full_path='I:\analysis\slices\2013-12-12-3-PH\data005'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

PrimaryNeurons=[ 257 289 588 616 679 1008 1097 1697 4293 4351     197  694 1235 1472  4399];
RecordingElectrodes=[ 10 15 40 42 46  68    66  114  291  291      22   47   83  107   294];
StimulatingElectrodeIDs=[ 4 1 1 1 3    3     3   1    4    4       4    1    2    2     4];
StimElectrodes=[127 103 174 199];

PrimaryNeurons=[ 257 289 588 616 679 1008 1097 1697 4293 4351     197  694 1235 1472  4399];
RecordingElectrodes=[ 10 15 40 42 46  68    66  114  291  291      22   47   83  107   294];
StimulatingElectrodeIDs=[ 4 1 1 1 3    3     3   1    4    4       4    1    2    2     4];



L=600;

for Neuron=1:length(PrimaryNeurons)
    figure(10)
    clf;
    
    fid=fopen(['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\dane_combined\ID=' num2str(PrimaryNeurons(Neuron)) 'c'],'r','ieee-le')
    a=fread(fid,'int32');
    l=length(a);
    b=reshape(a,5,l/5); %MovieNumber,RepetitionNumber,PatternNumber,Latency,TimeRelativeToRepetitionBegin
    fclose(fid);
    
    Electrode=StimulatingElectrodeIDs(Neuron)
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
            d1=d0(RecordingElectrodes(Neuron)+1,:);
            spikes(repetition,:)=d1;
        end
        subplot(4,2,movie);
        h=plot(spikes');
        set(h,'Color','b');
        hold on
        grid on
        
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
                d1=d0(RecordingElectrodes(Neuron)+1,:);
                h1=plot([LatencyOfSpike-10:LatencyOfSpike+20],d1([LatencyOfSpike-10:LatencyOfSpike+20]));
                %h1=plot(d1)
                set(h1,'Color','r');
            end
        end
        end                                                                
    end
    FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\report_2014_02_13\figures2\';
    FullImageName=[FigurePath 'Neuron' num2str(PrimaryNeurons(Neuron)) 'add.tif'];
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 9]);
    set(h,'PaperPosition',[0 0 16 9]); 
    print(h, '-dtiff', '-r120', FullImageName);
end    
