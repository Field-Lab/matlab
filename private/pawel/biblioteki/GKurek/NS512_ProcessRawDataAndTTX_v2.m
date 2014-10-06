function OutputData=NS512_ProcessRawDataAndTTX(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid)

%check source directory
if (full_path(length(full_path))=='\') separator=''; else separator='\'; end
files = dir([full_path separator 'data*.bin']);
if(size(files, 1)==0) error(['There are no files to process in ' full_path]); end

%check destination directory
if (output_path(length(output_path))=='\') separator=''; else separator='\'; end
files = dir([full_path separator 'data*.bin']);
if(size(files, 1)==0) error(['There are no files to process in ' output_path]); end

samplesInFile = 2400000;

if (full_path(length(full_path))=='\') shiftN=1; else shiftN=0; end
inputNumber = full_path(length(full_path)-2-shiftN:length(full_path)-shiftN);

fprintf('PROCESSING MOVIE NUMBER %d - [START]\n',MovieID);
fprintf(logFid,'PROCESSING MOVIE NUMBER %d - [START]\n',MovieID);

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);         

cd (full_path)
cd ..
NS_GlobalConstants = NS_GenerateGlobalConstants(512);
MovieData4=NS_MovieData(inputNumber,MovieID,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
fprintf(logFid,'-->Movie info:\n');
fprintf(logFid,'\tPDChunkNumber = %d\n', PDChunkNumber);
fprintf(logFid,'\tMovieBegin = %d\n', MovieBegin);
fprintf(logFid,'\tMovieEnd = %d\n', MovieBegin+ RepetNumber*RepetPeriod);
fprintf(logFid,'\tRepetNumber = %d\n', RepetNumber);
fprintf(logFid,'\tRepetPeriod = %d\n', RepetPeriod);
fprintf(logFid,'\n');
fprintf(logFid,'-->Processing details:\n');


fprintf('-->Movie info:\n')
fprintf('\tPDChunkNumber = %d\n', PDChunkNumber);
fprintf('\tMovieBegin = %d\n', MovieBegin);
fprintf('\tMovieEnd = %d\n', MovieBegin+ RepetNumber*RepetPeriod);
fprintf('\tRepetNumber = %d\n', RepetNumber);
fprintf('\tRepetPeriod = %d\n', RepetPeriod);
fprintf('\n');
fprintf('-->Processing details:\n')

ArtifactData = NS512_GetArtifactDataFromTTX(full_pathTTX, MovieID);

header = rawFile.getHeader();
binRep = header.getBinaryRepresentation();

files = dir([output_path separator 'data*.bin'])
fileName = files(1).name

NS_GlobalConstants = NS_GenerateGlobalConstants(512);
full_path

MovieData4=NS_MovieData_GlobalPath(MovieFilePath,MovieID,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
StimulusTimes=MovieData([1:3:length(MovieData)]);

for i=1:RepetNumber            
    startPosition=MovieBegin+(i-1)*RepetPeriod;
    FileNumber=floor(startPosition/samplesInFile);   
    
    RawData = rawFile.getData(startPosition,RepetPeriod)';    
     
    OutputData = RawData - ArtifactData; 
    for stimulus=1:length(StimulusTimes);
        OutputData(:,[StimulusTimes(stimulus)-1:StimulusTimes(stimulus)+16])=0;
    end

    out=NS512_CompressData16bit12bit(OutputData);    
    
    output_path_current=[output_path separator fileName(1:length(fileName)-6) num2str(FileNumber,'%02d') '.bin']
    fid = fopen(output_path_current, 'r+');
    offset = (startPosition - FileNumber*samplesInFile)*770;
    if FileNumber==0
        offset=offset+size(binRep,1);
    end
    fseek(fid, offset, 'bof');
    fwrite(fid, out, 'uint8');
    fclose(fid); 
    fprintf('\tRepetNumber:%d\t startPosition:%d\t FileNumber:%d\n', i, startPosition, FileNumber);
    fprintf(logFid,'\tRepetNumber:%d\t startPosition:%d\t FileNumber:%d\n', i, startPosition, FileNumber);

end
       fprintf('\n'); 
       fprintf(logFid,'\n'); 
       
if (MovieID ~= NumberOfMovies)
%if 2==5 %czy to ponizej naprawde jest potrzebne???
    
    MovieData4Next=NS_MovieData(inputNumber,MovieID+1,NS_GlobalConstants);
    [PDChunkNumberNext,MovieBeginNext,RepetNumberNext,RepetPeriodNext,MovieDataNext]=NS_DecodeMovieDataChunk(MovieData4Next);
   MovieDataNext
    startPosition = startPosition+RepetPeriod;
    fillWithZeros = MovieBeginNext - startPosition
    
    if (mod(fillWithZeros,10000)~=0)
        error('mod(fillWithZeros,10000)');
    end

    for i=1:(fillWithZeros/10000)
        
    startPosition=startPosition+(i-1)*RepetPeriod;
    FileNumber=floor(startPosition/samplesInFile);      
 
    OutputData = zeros(513,10000);      
   
    out=NS512_CompressData16bit12bit(OutputData);
    
    output_path_current=[output_path separator fileName(1:length(fileName)-6) num2str(FileNumber,'%02d') '.bin'];
    fid = fopen(output_path_current, 'r+');
    offset = (startPosition - FileNumber*samplesInFile)*770;
    if FileNumber==0
        offset=offset+size(binRep,1);
    end
    fseek(fid, offset, 'bof');
    fwrite(fid, out, 'uint8');
    fclose(fid);          
    fprintf('\tEraseDataAt:%d\t startPosition:%d\t FileNumber:%d\n', i, startPosition, FileNumber);
    fprintf(logFid,'\tEraseDataAt:%d\t startPosition:%d\t FileNumber:%d\n', i, startPosition, FileNumber);
    
    end

end

fprintf('PROCESSING MOVIE NUMBER %d - [END]\n',MovieID);
fprintf('\n');
fprintf('\n');

fprintf(logFid,'PROCESSING MOVIE NUMBER %d - [END]\n',MovieID);
fprintf(logFid,'\n');
fprintf(logFid,'\n');
end