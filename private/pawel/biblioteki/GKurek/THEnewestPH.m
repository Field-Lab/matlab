%clear
%clc
NS_GlobalConstants = NS_GenerateGlobalConstants(512);
cd G:\2010-09-14-0

fid0=fopen('movie002','r','b');
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies;
fclose(fid0);

full_path='G:\2010-09-14-0\data002\';           %INPUT
full_pathTTX='G:\2010-09-14-0\data009\';        %TTX INPUT
output_path = 'D:\Home\Pawel\analysis\2010-09-14-0\data002min009\data002000.bin';           %OUTPUT
%vision_path = 'E:\2010\Vision.jar';     %VISION

%javaaddpath(vision_path);   
samplesInFile = 2400000;

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);         
rawFileTTX=edu.ucsc.neurobiology.vision.io.RawDataFile(full_pathTTX); 

MovieID=13;

MovieData4=NS_MovieData('009',MovieID,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
PDChunkNumber;
MovieBegin;
RepetNumber;
RepetPeriod;

startPosition=MovieBegin;
Data = zeros(513,RepetPeriod,'int16'); %nazwa np. ArtifactData
for i=1:RepetNumber
    i
    RawDataTTX = rawFileTTX.getData(startPosition,RepetPeriod)';
    Data = int16(Data) + int16(RawDataTTX);
    startPosition=startPosition+RepetPeriod;
end

Data = Data/RepetNumber;

MovieData4=NS_MovieData('002',MovieID,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);

header = rawFile.getHeader();
binRep = header.getBinaryRepresentation();

startPosition=MovieBegin;
offset = size(binRep,1) + 770*MovieBegin;

%fid = fopen(output_path, 'r+');
for i=1:RepetNumber    
    startPosition=MovieBegin+(i-1)*RepetPeriod
    FileNumber=floor(startPosition/samplesInFile)    
    
    RawData = rawFile.getData(startPosition,RepetPeriod)';    
     
    OutputData = RawData - Data;        
    
    out=NS512_CompressData16bit12bit(OutputData);
    
    output_path_current=[output_path(1:length(output_path)-6) num2str(FileNumber,'%02d') '.bin']
    fid = fopen(output_path_current, 'r+');
    offset = (startPosition - FileNumber*samplesInFile)*770;
    if FileNumber==0
        offset=offset+size(binRep,1);
    end
    fseek(fid, offset, 'bof');
    fwrite(fid, out, 'uint8');
    fclose(fid);            
end

%fclose(fid);