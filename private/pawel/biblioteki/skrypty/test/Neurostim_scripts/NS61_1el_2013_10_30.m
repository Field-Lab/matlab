TimeShiftInMs=5;
DelayInMs=15;

NumberOfSamples=20000;
%electrodes=[1:8 10:24 26:56 58:64];
electrodes=[1:64];
Array=eye(64); 

electrodeOrder1 = LaurenelectrodeOrderGenerator(0);
electrodeOrder2 = LaurenelectrodeOrderGenerator(0);
Patterns=[electrodeOrder1];

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;
Times=[TimeShift:Delay:TimeShift+Delay*60];  % good only for 61 electrodes!!

Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie

cd C:\home\Pawel\nauka\NS61\stim_files; 
fid = fopen('NS61_1el_electrodes','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('NS61_1el_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('NS61_1el_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 