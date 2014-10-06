function [PatternOrder, MovieChunksFile]=NS_MovieChunksForExperimentSpatialPatterns(nElec, TimeShiftInMs, DelayInMs, nSamples)

% generates movieChunksFile for stimulus file generation
%
% specifies order in which patterns are played during stimulus (time each pattern is played)
%
% note: each movie chunk goes through all its repetitions (at one amplitude) before next movie chunk
% is played
%
% note: takes < 10s for stimulus to be generated when switching between movie chunks during
% experiment
%
% movie chunk length is limited to <=2 seconds, so patterns must be broken up into movie chunks, corresponding
% to different rows in Pattern
%
%%%%% info about corresponsing pattern array %%%%%%
%   pattern 1: all electrodes stimulated simultaneously at amplitudes that give 100% response rate
%   patterns 2-(1+nElec): each electrode stimulated individually at amplitude that gives 100% response
%   rate
%
%   patterns (2+nElec)-(1+6*nElec): all electrodes stimulated simultaneously, with all but one of the
%   electrodes stimulated at amplitudes that give 100% reponses, and last electrode stimulated at
%   amplitude at one of 5 amplitudes around threshold
%
%   patterns (2+6*nElec)-(1+11*nElec): each electrode stimulated individually at amplitude near
%   threshold
%
% order of patterns played:
% movie chunk 1:
%   first pattern: all electrodes simultaneously stimd at amp to get 100% response rate
%   patterns 2-(nElec*4+1): each electrode stimulated individually, at 1 ms offsets, in 4 different
%   orders
%   patterns (nElec*4+2)-(nElec*8+1): as previous but 5 ms offsets
%   patterns (nElec*8+2)-(nElec*12+1): as previous but 10 ms offsets
%
% movie chunks 2-nChunks: as movie chunk 1, but withe one of the electrodes at an amplitude near
% threshold instead of at 100% response level
%
% 2009-08-24, written by Lauren

randOrders = cell(4,1);
%randOrders{1} = randperm(nElec);
randOrders{1} = [1 2 3 4];

for i = 1:nElec %reverse order
    randOrders{3}(i) = randOrders{1}(nElec-i+1);
end
randOrders{2} = [randOrders{1}(2:end) randOrders{1}(1)];
randOrders{4} = [randOrders{3}(2:end) randOrders{3}(1)];

nChunks = 5*nElec+1;
nPatternsPerChunk = 12*nElec+1; %(4 orders)(3 delays)(nElecs) + simultaneous


Pattern = zeros(nChunks, nPatternsPerChunk);
% values in Pattern identify order in wich stimulus patterns (specified in patterns array, which is generated
% by a different function) are played



% pattern 1: all electrodes at amplitude necessary to elicit 100% response rate
% simultaneous
Pattern(1, 1) = 1;
% 1 ms delays
Pattern(1, 2:nElec*4+1) =          [1+randOrders{1}, 1+randOrders{2}, 1+randOrders{3}, 1+randOrders{4}];
% 5 ms delays
Pattern(1, nElec*4+2:nElec*8+1) =  [1+randOrders{1}, 1+randOrders{2}, 1+randOrders{3}, 1+randOrders{4}];
% 10 ms delays
Pattern(1, nElec*8+2:nElec*12+1) = [1+randOrders{1}, 1+randOrders{2}, 1+randOrders{3}, 1+randOrders{4}];


% all other patterns: all but one electrodes at amplitude necessary to elicit 100% response rate,
% remaining electrode at amplitude near threshold

offset = 6*nElec+1; %number of patterns before individual electrodes amps around threshold start

for i = 1:nElec
    for j = 1:5

        % simultaneous
        Pattern((i-1)*5+j+1, 1) = 1+nElec+5*(i-1)+j;

        patternsInMovieChunk = zeros(1, nElec); %pattern numbers corresponding to individual electrode
        % stimulation used in this movie chunk
        for k = 1:nElec
            if i == k
                patternsInMovieChunk(k) = offset+5*(i-1)+j;
            else
                patternsInMovieChunk(k) = k+1;
            end
        end

        % 1 ms delays
        Pattern((i-1)*5+j+1, 2:nElec*4+1) =          [patternsInMovieChunk(randOrders{1}) patternsInMovieChunk(randOrders{2})...
            patternsInMovieChunk(randOrders{3}) patternsInMovieChunk(randOrders{4})];

        % 5 ms delays
        Pattern((i-1)*5+j+1, nElec*4+2:nElec*8+1) =  [patternsInMovieChunk(randOrders{1}) patternsInMovieChunk(randOrders{2})...
            patternsInMovieChunk(randOrders{3}) patternsInMovieChunk(randOrders{4})];

        % 10 ms delays
        Pattern((i-1)*5+j+1, nElec*8+2:nElec*12+1) = [patternsInMovieChunk(randOrders{1}) patternsInMovieChunk(randOrders{2})...
            patternsInMovieChunk(randOrders{3}) patternsInMovieChunk(randOrders{4})];
    end
end

TimeShift = TimeShiftInMs*20; %in sampling periods (50 microseconds)
Delay = DelayInMs*20; %number of samples between patterns on different clusters

Times = zeros(1, nPatternsPerChunk);
Times(2:nElec:2+nElec*11) = Delay:Delay:Delay*12; %first time in each spatiotemporal pattern
Times(2:nElec:2+nElec*11) = Delay:Delay:Delay*3*nElec; %first time in each spatiotemporal pattern

for i = 1:nElec-1
    Times(2+i:nElec:2+i+nElec*3)          = Delay+20*i:Delay:Delay*4+20*i; %1 ms offsets
    Times(2+i+nElec*4:nElec:2+i+nElec*7)  = Delay*5+100*i:Delay:Delay*8+100*i; % 5 ms offsets
    Times(2+i+nElec*8:nElec:2+i+nElec*11) = Delay*9+200*i:Delay:Delay*12+200*i; % 10 ms offsets
end

Times = Times+TimeShift;

% makes a movie chunk for each row of Pattern and concatenates MovieChunks
MovieChunks = nChunks;
for i=1:nChunks
    Patterns = Pattern(i, :);
    Chunk = NS_MovieChunkGenerationForExperiment(Times, nSamples, Patterns);
    MovieChunks = [MovieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end

PatternOrder=Pattern;
MovieChunksFile=MovieChunks;

