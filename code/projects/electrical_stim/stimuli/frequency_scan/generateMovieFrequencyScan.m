function MovieChunks=generateMovieFrequencyScan(PatternNumber, Periods, DelayInMs, NumberOfSamples, maxReps)

% generates a stimulus movie that includes as many pulses as can fit into a single movie chunk at
% each frequence, but limits maximum number of pulse applications at a given frequency to the value
% maxReps
%
% patternNumber: number corresponding to pattern that is applied at different frequencies
% (typically, this is just the number of the electrode used for stimulation)
%
% periods: set of periods (in samples) corresponding to frequencies tested, in the order they will
% be applied
%
% DelayInMs: time delay before movie chunk is started after each trigger (typically 0--code untested
% for nonzero values)
%
% NumberOfSamples: length of time for each movie chunk, in samples (can't be more than 40000)
%
% maxReps: number of stimulus applications for each frequency
% 
% code checked October 2010

delay = DelayInMs*20;

MovieChunks=length(Periods); %one movie chunk for each frequency
for i=1:length(Periods)
    if (NumberOfSamples-delay-100)/Periods(i) > maxReps
        Times=delay:Periods(i):(delay+Periods(i)*(maxReps-1)); %do maxReps repetitions at this frequency
    else
        Times=delay:Periods(i):(NumberOfSamples-Periods(i)); %do as many repetitions at this frequency as fit into movie chunk length
    end
    
    Patterns=ones(1,length(Times))*PatternNumber;
    
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
    MovieChunks=[MovieChunks Chunk]; %#ok<AGROW>
end