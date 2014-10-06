clear
clc
NS_GlobalConstants = NS_GenerateGlobalConstants(512);
cd J:\2010-09-14-0
full_path='J:\2010-09-14-0\data002\';           %INPUT
output_path = 'J:\2010-09-14-0\';       %OUTPUT
vision_path = 'E:\2010\Vision.jar';     %VISION

javaaddpath(vision_path);  

MovieID=6;
MovieData4=NS_MovieData('002',MovieID,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);


rawFileIN=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);   
rawFileOUT=edu.ucsc.neurobiology.vision.io.RawDataFile(output_path);   

header = rawFileIN.getHeader();
binRep = header.getBinaryRepresentation();

RawDataIN= rawFileIN.getData(2800000,60000)';
RawDataOUT= rawFileOUT.getData(2800000,60000)';
clf
plot(RawDataIN(100,:));
hold on
plot(RawDataOUT(100,:),'r');