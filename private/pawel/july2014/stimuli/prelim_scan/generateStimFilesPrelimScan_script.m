% note: on 2010-03-02 this script was changed from playing only one random order of electrodes to
% playing 2 random orders of electrodes

%SET TRIGGER INTERVAL TO 1 SECOND

cd('/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test')

TimeShiftInMs=0;
DelayInMs=7.5; %time between subsequent pulses on different electrodes

NumberOfSamples = 20000;
% refresh rate (trigger interval) in labview must be the same as this number!!!
electrodes=1:64;
Array=eye(64);

electrodeOrder1 = LaurenelectrodeOrderGenerator(0);
electrodeOrder2 = LaurenelectrodeOrderGenerator(0);
Patterns=[electrodeOrder1 electrodeOrder2];

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;
Times=(0:Delay:Delay*121) + TimeShift; %121 = length(Patterns) - 1

Chunk=NS_MovieChunkGenerationForExperiment(Times, NumberOfSamples, Patterns);
MovieChunksFile=[1 Chunk]; %only one movie


fid = fopen('1el_prel_electrodes','wb','ieee-le.l64')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('1el_prel_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('1el_prel_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);