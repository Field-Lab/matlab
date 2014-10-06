% code checked October 2010

cd '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_experiment'

%% SET TRIGGER INTERVAL TO 1 s

PatternNumber=53; %electrode to apply frequency scan on




TimeShiftInMs=0;
DelayInMs=0;
NumberOfSamples=20000;

Periods=[200 100 50 25 12.5 6.25]*20; %in samples
maxReps = 20; %limits number of pulses at each frequency (for a single labview iteration)

electrodes=1:64;
Array=eye(64);

for jj = 1:length(Periods)
    
    period = Periods(jj);
    MovieChunksFile=generateMovieFrequencyScan(PatternNumber, period, DelayInMs, NumberOfSamples, maxReps);

    electrodes=1:64;
    Array=eye(64);

    
    prefix = ['freq_scan_' num2str(PatternNumber) '_freq_' num2str(20000/period)];

    fid = fopen([prefix '_electrodes'],'wb','ieee-le.l64')
    fwrite(fid,electrodes,'int32');
    fclose(fid);

    fid = fopen([prefix '_patterns'],'wb','ieee-le.l64')
    fwrite(fid,Array,'double');
    fclose(fid);

    fid = fopen([prefix '_movie'],'wb','ieee-le.l64')
    fwrite(fid,MovieChunksFile,'int32');
    fclose(fid);
end



