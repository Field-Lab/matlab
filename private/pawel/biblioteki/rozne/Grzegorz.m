clear
clc

%-Paths---
full_path='D:\2010\data003\';     %define path to raw data file
full_pathTTX='D:\2010\data005\';  %define path to raw data file (TTX)
output_path = 'D:\2010\out\';       %define path to save results
    

%-Variables---
javaaddpath 'E:\2010\Vision.jar';               %define path to Vision jar file
samplesInFile = 1200000;

files    = length(dir([full_path 'data*.bin']));
filesTTX = length(dir([full_pathTTX 'data*.bin']));

if (files==0) 
    error( [full_path '  <--There is no data to process.']); 
end;

if (filesTTX==0) 
    error([full_pathTTX '  <--There is no TTX data to process.']); 
end;

amountOfFiles = min(files,filesTTX);
currentFile = 0;
dataToProcess = samplesInFile * amountOfFiles;
output_file = [output_path 'data00' '_' int2str(currentFile) '.bin'];



rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);         
rawFileTTX=edu.ucsc.neurobiology.vision.io.RawDataFile(full_pathTTX); 

header = rawFile.getHeader();
binRep = header.getBinaryRepresentation();


%-Save header
fid = fopen(output_file, 'w');
if (fid<0)
    error('Check output dirrectory');
end

fwrite(fid, binRep, 'int8');
    
StartPos = 0;           %Start position
SampleAm = 2000;        %Amount of data

tic

fprintf('0.00%% done\n');
fprintf('Elapsed time is 00.000000 seconds.');

for k=StartPos:SampleAm:(dataToProcess-SampleAm),
    
    if( (k~=0) && (~mod(k,samplesInFile)) )
        fclose(fid);
        output_file = [output_path 'data00' '_' int2str(k/samplesInFile) '.bin'];
        fid = fopen(output_file, 'w');
    end
    
    %-Get specifed data---
    RawData=rawFile.getData(k,SampleAm)';   
    RawDataTTX=rawFileTTX.getData(k,SampleAm)';
    
    %-Result-difference between RawData and Data with TTX---
    OutputData = RawData - RawDataTTX; 
    

    for i=1:SampleAm,
        fwrite(fid, OutputData(1,i), 'int16');
        vector = [0,0];
        position=1;
        
        for j=2:2:512,
            
            s1 = uint16( OutputData(j,i)+2048 );
            s2 = uint16( OutputData(j+1,i)+2048 );
            
            
            b1 = bitshift(s1,-4);
            %fwrite(fid, b1, 'uint8');
            vector(position)=b1;
            position=position+1;
            
            b2 = bitshift(bitand(s1,15),4) + bitshift(s2,-8);
            %fwrite(fid, b2, 'uint8');
            vector(position)=b2;
            position=position+1;
            
            b3=bitand(s2,255);
            %fwrite(fid, b3, 'uint8');
            vector(position)=b3;
            position=position+1;
  
        end
        
        fwrite(fid, vector, 'uint8');

        
    end

           clc
        fprintf('%1.2f%% done\n',k/(dataToProcess-SampleAm)*100);
        toc
end

fclose(fid);


    
    
