TimeShiftInMs=10;
DelayInMs=100;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512); 

Patterns=[];
Times=[];
Electrode=25; %change this!

Chunk=NS_MovieChunkGenerationForExperiment(2000,NumberOfSamples,Electrode);
MovieChunksFile=[1 Chunk]; 
break
cd C:\home\pawel\nauka\SantaCruz\Stim_system_testing; %change this path
fid = fopen('Few_el_electrodes2','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('Few_el_patterns2','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('Few_el_movie2','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 