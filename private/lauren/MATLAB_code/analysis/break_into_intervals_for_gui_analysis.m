function writtenData = break_into_intervals_for_gui_analysis(pathToDataFile, WritePath, interval, movieRepsPerMovieOut, movieNoTimes, patternNo, patternStrOrig, origAnalysisPath)


% splits an extended recording period into intervals that can be analyzed
% using the standard analysis gui
% interval is the length between break points, but an extra 10 samples is
% saved in overlap


%ArrayID=1;

if ~exist([WritePath filesep 'pattern_files'], 'file')
    mkdir([WritePath filesep 'pattern_files'])
end

if ~exist([WritePath filesep 'status_files'], 'file')
    mkdir([WritePath filesep 'status_files'])
end


rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(pathToDataFile);

totalSamples = getRawDataFileSize(pathToDataFile);

writtenData = false(1,totalSamples);


nMoviesOut = size(movieNoTimes,2)/movieRepsPerMovieOut;

%number of times the original movie is applied must be divisible by
%movieRepsPerMovieOut
nRepsOriginal = sum(movieNoTimes(1,:)==movieNoTimes(1,1));
if mod(nRepsOriginal,movieRepsPerMovieOut) ~= 0
    error('original number of repetitions of pattern at each amplitude must be divisible by movieRepsPerMovieOut')
end    

nRepsPerMovieRepIn = (movieNoTimes(2,2)-movieNoTimes(2,1))/interval;
%nMovies = size(movieNoTimes,2)*(movieNoTimes(2,2)-movieNoTimes(2,1))/(intsPerMovie*interval);

Channels=1:64;
traceLength = interval + 10;

nReps = nRepsPerMovieRepIn*movieRepsPerMovieOut; %total number of reps to be packaged into each movie out
repOutStartTimes = zeros(nMoviesOut, nReps); %sample time, relative to pulse applied in original movie, corresponding to first sample of repetition in output file

% iterates through movies (movies chunks x stimulus amplitudes)
for ii = 1:nMoviesOut
    disp(ii)
    
    Traces = int16(zeros(nReps, length(Channels), traceLength)); %format of tracesToSave = reps x channels x samples
        
    for jj = 1:movieRepsPerMovieOut 
        movieRepInd = (ii-1)*movieRepsPerMovieOut + jj;
        origMovNo = movieNoTimes(1,movieRepInd);
        startSample = movieNoTimes(2,movieRepInd);
        for kk = 1:nRepsPerMovieRepIn;
            if ~(ii == nMoviesOut && jj == movieRepsPerMovieOut && kk == nRepsPerMovieRepIn)
                RawData=int16(rawFile.getData(startSample+(kk-1)*interval, traceLength)');
                writtenData(startSample+(kk-1)*interval:startSample+(kk-1)*interval+traceLength) = true;
            else
                RawData = zeros(65,traceLength);
                RawData(:,1:interval-1) = int16(rawFile.getData(startSample+(kk-1)*interval, interval-1)');
                writtenData(startSample+(kk-1)*interval:startSample+(kk-1)*interval+interval) = true;
            end
            Traces((jj-1)*nRepsPerMovieRepIn + kk, :, :) = RawData(Channels+1,:);
        
            repOutStartTimes(ii, (jj-1)*nRepsPerMovieRepIn + kk) = (kk-1)*interval;
        end
    end
        
    %save the pattern/movie file
    STTS=size(Traces);
    
    a=reshape(Traces,STTS(1)*STTS(2)*STTS(3),1);
    b=zeros(1000,1);
    b(1)=STTS(1);
    b(2)=STTS(2);
    b(3)=STTS(3);
    b(3+1:3+length(Channels))=Channels';
    o=[b' a'];
    
    pString = ['p' num2str(patternNo) '_om' num2str(origMovNo)];
    
    if ~exist([WritePath filesep pString], 'file')
        mkdir([WritePath filesep pString])
    end
    
    FullName=[WritePath filesep pString filesep pString '_m' num2str(ii)];
    fid=fopen(FullName,'wb','ieee-le');
    fwrite(fid,o,'int16');
    fclose(fid);
       
    %copy original pattern and status files to output folder under new corresponding names
    originalPatternFile = [origAnalysisPath filesep 'pattern_files' filesep 'pattern' patternStrOrig '_m' num2str(origMovNo) '.mat'];
    newPatternFile = [WritePath filesep 'pattern_files' filesep 'pattern' num2str(patternNo) '_om' num2str(origMovNo) '_m' num2str(ii) '.mat'];
    copyfile(originalPatternFile, newPatternFile)
    
    originalStatusFile = [origAnalysisPath filesep 'status_files' filesep 'status_m' num2str(origMovNo) '.mat'];
    newStatusFile = [WritePath filesep 'status_files' filesep 'status_m' num2str(ii) '.mat'];
    copyfile(originalStatusFile, newStatusFile)
end
    

save([WritePath filesep 'repStartTimes.mat'], 'repOutStartTimes')

end

