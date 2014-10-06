function y=NS_ShowChannelDataExtended(ChannelsInfo,MoviesInfo,FileName,Channel,AdditionalChannels,TimeStart,TimeLength,Unit,FigureNumber,NS_GlobalConstants)
%A simple function that shows the raw data and the stimulation current
%waveforms on the saeme figure.
%FileName - format like: '003';

Fs=NS_GlobalConstants.SamplingFrequency;
time=[TimeStart+1 TimeStart+TimeLength];
t=[TimeStart+1:TimeStart+TimeLength]/Fs;
%Find the movies that include the data of interest
movies=[];
for i=1:length(MoviesInfo)
    start=MoviesInfo(i).StartTime;
    movie_length=MoviesInfo(i).PeriodOfRep*MoviesInfo(i).NumberOfRep;
    if (start+1<time(2) & start+movie_length>time(1))
        movies=[movies i];
    end
end
movies
channels=[Channel];
ChannelsInfo_size=size(ChannelsInfo);
for i=1:length(movies)
    pattern=ChannelsInfo(Channel,movies(i));
    if pattern~=0
        for j=1:ChannelsInfo_size(1)
            if ChannelsInfo(j,movies(i))==pattern
                if (length(find(channels==j))==0 & j~= Channel)
                channels=[channels j];
                end
            end
        end
    end
end

channels=[channels AdditionalChannels]

filename_movie=['movie' FileName];

figure(FigureNumber);
for i=1:length(channels)
    output1=ReadMovieDataChunk4(filename_movie,channels(i),time,Unit,NS_GlobalConstants);
    subplot(6,length(channels),i);
    plotPH(t,output1(2,:),100);
    subplot(6,length(channels),length(channels)+i);
    plotPH(t,output1(4,:),100);
    subplot(3,length(channels),length(channels)+i);
    plotPH(t,output1(1,:).*output1(3,:),100);
    clear output1;
    y=gca;
    xlabel('t');
    %ylabel('I')
end

full_path=[pwd '\' 'data' FileName];
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
data = rawFile.getData(TimeStart+1,TimeLength);
for i=1:length(channels)
    subplot(3,length(channels),2*length(channels)+i);
    plot(t,data(:,channels(i)+1));
    y=gca;
    xlabel('t')
    %ylabel('I')
end

y=channels;