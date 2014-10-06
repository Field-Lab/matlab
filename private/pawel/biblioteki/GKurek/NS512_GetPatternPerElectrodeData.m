function [electrodeRawData, patternRepeatData] = NS512_GetPatternPerElectrodeData(movieID, patternID, electrodeID, full_path)

currentWorkspace = pwd;

cd (full_path)
cd ..
NS_GlobalConstants = NS_GenerateGlobalConstants(512);
MovieData4=NS_MovieData('002',movieID,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
cd (currentWorkspace)

patternIDList = MovieData(2:3:end,:);
patternTimeList = MovieData(1:3:end,:);

patternIDPositions = find(patternIDList==patternID);

if(size(patternIDPositions,1)==0) 
    fprintf('Patern number=%d not found in movie number=%d\n', patternID, movieID); error('-1');
end

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);         

patternExecutionTime = 600;

electrodeRawData = [];
patternRepeatData = [];

for i=1:size(patternIDPositions,1)  
    
    patternIDStartTimeInMovie = patternTimeList(patternIDPositions(i));
    patternIDStartTimeGlobal = patternIDStartTimeInMovie + MovieBegin;
    %figure;
    hold on;
    for j=0:RepetNumber-1
        RawData = rawFile.getData(patternIDStartTimeGlobal+j*RepetPeriod,patternExecutionTime)';
        electrodeRawData = [electrodeRawData ; RawData(electrodeID,:)];
        patternRepeatData = [patternRepeatData;i];
    end    
end

end

