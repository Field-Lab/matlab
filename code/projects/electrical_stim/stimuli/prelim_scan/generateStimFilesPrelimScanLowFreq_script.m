% note: on 2010-03-02 this script was changed from playing only one random order of electrodes to
% playing 2 random orders of electrodes

%SET TRIGGER INTERVAL TO 1.6 SECONDS
%4 MOVIE CHUNKS
%ONLY 1 RANDOM ORDER SO SET REPS IN LABVIEW TO 52!!!
%
%requires ~5 minutes/amplitude so choose small range or bigger steps (for
%15 amplitudes with 1.1x steps, upper end of range should be 4x lower end)

cd('/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test')

TimeShiftInMs=0;
DelayInMs=100; %time between subsequent pulses on different electrodes

%NumberOfSamples = DelayInMs*20*61*2; %(delay between pulses)(conversion to samples)(61 stimuli)(2 orders)
NumberOfSamples = 32000;
% refresh rate (trigger interval) in labview must be the same as this number!!!
electrodes=1:64;
Array=eye(64);

%only use one random order because they need to be broken into multiple
%movie chunks
Patterns = LaurenelectrodeOrderGenerator(0);

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;
%Times=(0:Delay:Delay*60) + TimeShift; %121 = length(Patterns) - 1



%break into 1-16, 17-32, 33-48, 49-61

Pattern_chunks{1} = Patterns(1:16);
Pattern_chunks{2} = Patterns(17:32);
Pattern_chunks{3} = Patterns(33:48);
Pattern_chunks{4} = Patterns(49:61);

MovieChunksFile = 4; %4 chunks
for ii = 1:4
    Times = 0:Delay:Delay*(length(Pattern_chunks{ii})-1) + TimeShift;
    Chunk = NS_MovieChunkGenerationForExperiment(Times, NumberOfSamples, Pattern_chunks{ii});
    MovieChunksFile = [MovieChunksFile Chunk];
end

keyboard


fid = fopen('1el_low_freq_electrodes','wb','ieee-le.l64')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('1el_low_freq_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('1el_low_freq_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);