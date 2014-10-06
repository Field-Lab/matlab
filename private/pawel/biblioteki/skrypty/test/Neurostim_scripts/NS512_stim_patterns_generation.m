TimeShiftInMs=10;
DelayInMs=50;

NumberOfSamples=10000;
%electrodes=[1:8 10:24 26:56 58:64];
electrodes=[1:512];
%Array=eye(512); 

electrodeOrder1 = LaurenelectrodeOrderGenerator(0);
electrodeOrder2 = LaurenelectrodeOrderGenerator(0);
Patterns=[1:512];

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;
Times=[TimeShift:Delay:TimeShift+Delay*(length(Patterns)-1)];  % good only for 61 electrodes!!

%[electrodes,Array,Times]=NS512_LongMovingBarStimulus(10,30);

Columns=[63];
[electrodes,Array,Times]=NS512_ShortMovingBarStimulus(10,0,Columns);

Times=[TimeShift:Delay:NumberOfSamples];

Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,[1 1 1 1 1 1 1 1 1 1]);
MovieChunksFile=[1 Chunk]; %only one movie

%break;
cd C:\home\pawel\praca\stim_files; 
fid = fopen('NS512_ShortBar_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('NS512_ShortBar_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('NS512_ShortBar_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 