MovieFilePath='D:\Home\Data\slices\2010-09-14-0\movie002';
movieIDList = NS512_GetMoviesWithPatternIDNew(MovieFilePath,patternID);
MovieData4=NS_MovieData_GlobalPath(MovieFilePath,movieIDList(WhichMovie),NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
SMD=size(MovieData);
PatternNumbers=MovieData(2:3:SMD(1));
WhichPatterns=find(PatternNumbers==patternID);
WaveformLength=600;
spikes3=zeros(RepetNumber,WaveformLength);
figure(4)
for i=1:length(WhichPatterns)
    TimeIndex=MovieData(1+(WhichPatterns(i)-1)*3);
    PulseTimes=MovieBegin+TimeIndex+(0:RepetNumber-1)*RepetPeriod;
    
    for t=1:RepetNumber
        PulseTime=PulseTimes(t);
        d0=rawFile.getData(PulseTime,WaveformLength)';
        d1=d0(PrimaryElectrode+1,:);
        size(d1);
        spikes3(t,:)=d1;
        %spikes2(i,:,:)=d1;
    end
    subplot(1,length(WhichPatterns),i);
    plot(spikes3');
end