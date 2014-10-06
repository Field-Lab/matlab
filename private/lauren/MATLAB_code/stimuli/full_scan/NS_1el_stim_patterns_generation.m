%% SET TRIGGER INTERVAL TO 2 SECONDS

% cd('/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test')

TimeShiftInMs=0;
DelayInMs=15;

NumberOfSamples=40000;
electrodes=[1:64];
Array=eye(64); 

electrodeOrder1 = LaurenelectrodeOrderGenerator(0);
electrodeOrder2 = LaurenelectrodeOrderGenerator(0);
Patterns=[electrodeOrder1 electrodeOrder2];

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;


Times = (0:Delay:Delay*(length(Patterns)-1)) + TimeShift;
Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie

% cd C:\pawel\pliki\nauka\matlab\; 
%cd /Applications/MATLAB74/work/Lauren/stimuli/stimulus_files/; 
keyboard; 

fid = fopen('1el_electrodes','wb','ieee-le.l64')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('1el_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('1el_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 