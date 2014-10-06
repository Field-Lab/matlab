function [ArtifactData] = NS512_GetArtifactDataFromTTX( full_pathTTX, MovieID );

%check source directory
if (full_pathTTX(length(full_pathTTX))=='\') separator=''; else separator='\'; end
files = dir([full_pathTTX separator 'data*.bin']);
if(size(files, 1)==0) display(['There are no files to process in ' full_pathTTX]); end

if (full_pathTTX(length(full_pathTTX))=='\') shiftNTTX=1; else shiftNTTX=0; end
inputNumberTTX = full_pathTTX(length(full_pathTTX)-2-shiftNTTX:length(full_pathTTX)-shiftNTTX);

rawFileTTX=edu.ucsc.neurobiology.vision.io.RawDataFile(full_pathTTX);

cd (full_pathTTX)
cd ..
NS_GlobalConstants = NS_GenerateGlobalConstants(512);
MovieData4=NS_MovieData(inputNumberTTX,MovieID,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);

startPosition=MovieBegin;
Data = zeros(513,RepetPeriod,'int32'); %nazwa np. ArtifactData
for i=1:RepetNumber    
    RawDataTTX = int32(rawFileTTX.getData(startPosition,RepetPeriod))';
    Data = Data + RawDataTTX;
    startPosition=startPosition+RepetPeriod;
end

ArtifactData = int16(Data/RepetNumber);

end