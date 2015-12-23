function fakePreprocessStimData(FileName, WritePath, traceLength, timeRegionMs)

%makes fake partial status file
Status.ChannelsStatus = struct('range', ones(1,64));
FullName_status=[WritePath filesep 'status_m' num2str(999)];
save(FullName_status, 'Status')


full_path=[pwd filesep 'data' FileName] %#ok<NOPRT>
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

totalSamples = getRawDataFileSize(full_path);


if totalSamples <= timeRegionMs(2)*20;
    warndlg(['Raw data file ' full_path ' does not have as many samples as expected by stimulus files']);
    return
end

timeLength = timeRegionMs(2) - timeRegionMs(1) + 1;
nTraces = floor(timeLength*20/traceLength);


%filename_movie=['movie' FileName] %#ok<NOPRT>
%l=length(filename_movie);
%SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)] %#ok<NOPRT>

%fid0=fopen(filename_movie,'r','b');
%header=readMHchunk(fid0);
%nMovies=header.number_of_movies;

Channels=1:64;

TracesToSave=zeros(nTraces, length(Channels), traceLength, 'int16');

RawData=int16(rawFile.getData(timeRegionMs(1)*20, timeLength*20)');

for j=1:nTraces %iterates through number of times movie is played
    t = (timeRegionMs(1)-1)*20 + (j-1)*traceLength;
    TracesToSave(j,:,:)=RawData(Channels+1, t+1:t+traceLength);
end


STTS=size(TracesToSave);

a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
b=zeros(1000,1);
b(1)=STTS(1);
b(2)=STTS(2);
b(3)=STTS(3);
b(3+1:3+length(Channels))=Channels';
o=[b' a'];


FullName=[WritePath filesep 'p999_m999'];
fid=fopen(FullName,'wb','ieee-le');
fwrite(fid,o,'int16');
fclose(fid);

% FullName_pattern=[WritePath filesep 'pattern_files' filesep 'pattern' subgroupInfo(j).patternName(2:end) '_m' num2str(i)];
% save(FullName_pattern,'pattern');


clear TracesToSave;