TimeShiftInMs=10;
DelayInMs=100;

NumberOfSamples=10000;
%electrodes=[1:8 10:24 26:56 58:64];
electrodes=[1:512];
Array=eye(512); 

Patterns=[];
Times=[];
chip=3;
Electrodes=[1:64]+(chip-1)*64;
Electrodes=[1:8:512];

for i=1
    for el=1:length(Electrodes)
        Patterns=[Patterns Electrodes(el)];
        Times=[150:150:64*150];
    end        
end

Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie
break
cd C:\home\Pawel\nauka\StimFiles;
fid = fopen('1el_scan_electrodes','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('1el_scan_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('1el_scan_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 