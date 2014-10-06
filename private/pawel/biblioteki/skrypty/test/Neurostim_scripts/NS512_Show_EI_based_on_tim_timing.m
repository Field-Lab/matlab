NS_GlobalConstants=NS_GenerateGlobalConstants(512);

FolderName='G:\2010-09-11-0';
DataFilename='003';

full_path='G:\2010-09-11-0\data003';
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

cd(FolderName);

Info=NS512_MovieNumberAndTimingForPattern('G:\2010-09-11-0\movie003',508);
MovieNumber=Info(1,1);

XStart=Info(5,3);
RepetPeriod=10000;
Channels=[1:512]; %[1:64 321: 512];
Length=300;
Data=int16(zeros(length(Channels),Length));

ChunkData=NS_MovieData('003',MovieNumber,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);

for j=1:50
        j
        RawData=int16(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod)');  
        a=RawData(Channels+1,XStart:XStart+Length-1);
        Data(:,:)=Data(:,:)+a;           
end

Data=Data/50;