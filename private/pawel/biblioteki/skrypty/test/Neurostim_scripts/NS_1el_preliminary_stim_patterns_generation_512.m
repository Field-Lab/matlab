TimeShiftInMs=0;
DelayInMs=3;

NumberOfSamples=10000;
%electrodes=[1:8 10:24 26:56 58:64];
electrodes=[1:128]+128;
Array=eye(128); 

electrodeOrder1 = LaurenelectrodeOrderGenerator(0);
electrodeOrder2 = LaurenelectrodeOrderGenerator(0);
Patterns=[1:128];

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;
Times=[TimeShift:Delay:TimeShift+Delay*(length(Patterns)-1)];  % good only for 61 electrodes!!

Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie

cd C:\home\pawel\nauka\Neurostim512\scripts; 
fid = fopen('1el_prel_electrodes','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('1el_prel_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('1el_prel_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 