% note: on 2010-03-02 this script was changed from playing only one random order of electrodes to
% playing 2 random orders of electrodes

TimeShiftInMs=0;
DelayInMs=15; %time between subsequent pulses on different electrodes

NumberOfSamples = DelayInMs*20*61*2; %(delay between pulses)(conversion to samples)(61 stimuli)(2 orders)
% refresh rate in labview must be the same or larger than this number!!!
electrodes=1:64;
Array=eye(64);

electrodeOrder1 = LaurenelectrodeOrderGenerator(0);
electrodeOrder2 = LaurenelectrodeOrderGenerator(0);
Patterns=electrodeOrder1;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;
Times=(0:Delay:Delay*60) + TimeShift;


Chunk=NS_MovieChunkGenerationForExperiment(Times, NumberOfSamples, Patterns);
MovieChunksFile=[1 Chunk]; %only one movie


fid = fopen('1el_prel_old_electrodes','wb','ieee-le.l64')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('1el_prel_old_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('1el_prel_old_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 