function MovieChunksFile = generateMovieSpatiotempProbeExtended(periodInMs, offsets, details)

% generates movieChunksFile for stimulus file generation
%
% offsets: a vector of offset values specifying delays between subsequent pulses within a sequence
% (in ms)
% specifies order in which patterns are played during stimulus (time each pattern is played)
%
% movie chunk length is limited to <=2 seconds (nSamples <= 40000)
%
% generates 2 random orders of stimuli
%
%
% 2010-03-19: added second repetition of center-alone stimulus in final chunk so that total reps
% match other stimuli (since other stimuli have 2x random orders)
%
% note: in order for preprocessing to work correctly (correct parsing of
% raw data), each center electrode pulse amplitude must be in a separate
% movie chunk



prepulseBin = ~~details.prePulse;
stimPulseAmps = details.stimPulseAmps;
pSimAmps = details.pairsCenterAmps;
%clusterID = details.clusterID;

nOffsets = length(offsets);
period = periodInMs*20;

nPrepulseAmpCombos = sum(prepulseBin); %number of prepulse/pre-amplitude combinations
nSimPulseCombosPerPrimAmp = sum(pSimAmps == 1);


samplesPerPrimAmpOffset = nPrepulseAmpCombos*2*period; %factor of 2 for 2 random orders

chunksPerOffset = ceil(samplesPerPrimAmpOffset/40000); %chunks required to apply all prepulse-stimpulse combinations for a given stimpulse amp and offset

nPulsePairsPerChunk = ceil(nPrepulseAmpCombos/chunksPerOffset);

if nPulsePairsPerChunk*2*period+100 > 40000 %may happen if 40000 doesn't divide into periods well???
    keyboard
    chunksPerOffset = chunksPerOffset+1;
    nPulsePairsPerChunk = ceil(nPrepulseAmpCombos/chunksPerOffset);
end

nSamples = ceil((100+nPulsePairsPerChunk*2*period)/2000)*2000; %rounds up to the nearest 0.1 seconds, 100 added so that there is at least 100 samples after final pulse in chunk

warndlg(['set trigger interval to ' num2str(nSamples/20000) ' seconds'])


%chunks required to apply all simultaneous stim + 1 isolated stim for a given stimpulse amp
chunksPerPrimAmpSimPulses = ceil((period*2*(nSimPulseCombosPerPrimAmp + 1)/nSamples)); %factor of 2 for 2 random orders, +1 for isolated stim
nSimPulsesPerChunk = ceil(nSimPulseCombosPerPrimAmp/chunksPerPrimAmpSimPulses);

%total number of chunks required
nChunks = (nOffsets*chunksPerOffset + chunksPerPrimAmpSimPulses)*max(stimPulseAmps);

warndlg([num2str(nChunks) ' movie chunks required'])


prepulseInd = find(prepulseBin);
prepulseIndRand = prepulseInd(randperm(length(prepulseInd)));
%randomly pick with prepulses are divided into which chunks
for ii = 1:chunksPerOffset
    if ii*nPulsePairsPerChunk <= nPrepulseAmpCombos
        prepulseIndChunks{ii} = prepulseIndRand((ii-1)*nPulsePairsPerChunk+1:ii*nPulsePairsPerChunk);
    else
        prepulseIndChunks{ii} = prepulseIndRand((ii-1)*nPulsePairsPerChunk+1:end);
    end
end

simPulseOrderRand = randperm(nSimPulseCombosPerPrimAmp);
for ii = 1:chunksPerPrimAmpSimPulses
    if ii*nSimPulsesPerChunk <= nSimPulseCombosPerPrimAmp
        simOrderChunks{ii} = simPulseOrderRand((ii-1)*nSimPulsesPerChunk+1:ii*nSimPulsesPerChunk);
    else
        simOrderChunks{ii} = simPulseOrderRand((ii-1)*nSimPulsesPerChunk+1:end);
    end
end


iChunk = 0;
for ii = 1:max(stimPulseAmps)
    centerAloneInd = find(stimPulseAmps == ii);
    
    
    for jj = 1:nOffsets
        for kk = 1:chunksPerOffset
            
            iChunk = iChunk + 1;

            timesTemp = period:period:(period + period*2*length(prepulseIndChunks{kk}) - 1); %check
            
            %prepulses
            pattern{iChunk} = [prepulseIndChunks{kk}(randperm(length(prepulseIndChunks{kk})))...
                prepulseIndChunks{kk}(randperm(length(prepulseIndChunks{kk})))]; %2 random orders of prepulses
            times{iChunk} = timesTemp - offsets(jj)*20;
            
            %stim pulses
            pattern{iChunk} = [pattern{iChunk} centerAloneInd*ones(1,2*length(prepulseIndChunks{kk}))]; %stim pulses
            times{iChunk} = [times{iChunk} timesTemp];
        end
    end
    
    %last chunks for ith stimPulseAmp = simultaneous-pulse stimuli and the
    %center-alone stimuli
    
    simInd = find(pSimAmps == ii);
    
    for jj = 1:chunksPerPrimAmpSimPulses
        iChunk = iChunk+1;
        if jj < chunksPerPrimAmpSimPulses
            pattern{iChunk} = [simInd(simOrderChunks{jj}(randperm(length(simOrderChunks{jj}))))...
                simInd(simOrderChunks{jj}(randperm(length(simOrderChunks{jj}))))]; %#ok<AGROW> %first random order
            times{iChunk} = period:period:period*2*(length(simOrderChunks{jj}));
        else %add isolated stim pulse to final movie chunk
            pattern{iChunk} = [centerAloneInd simInd(simOrderChunks{jj}(randperm(length(simOrderChunks{jj}))))...
                centerAloneInd simInd(simOrderChunks{jj}(randperm(length(simOrderChunks{jj}))))]; %#ok<AGROW> %first random order
            times{iChunk} = period:period:period*2*(length(simOrderChunks{jj})+1);
        end
    end
end

if length(pattern) ~= nChunks
    error('number of chunks used does not match number of estimated required number of chunks, aborting!')
end



%% checks that each pattern is applied the expected number of times

patternsPlayed = zeros(1, length(stimPulseAmps));
for ii = 1:nChunks
    for jj = 1:length(pattern{ii})
        patternsPlayed(pattern{ii}(jj)) = patternsPlayed(pattern{ii}(jj)) + 1;
    end
end

for jj = 1:length(patternsPlayed)
    if prepulseBin(jj) && patternsPlayed(jj) ~= 2*nOffsets*max(stimPulseAmps)
        error(['unexpected number of applications of pattern ' num2str(jj) '... aborting'])
    elseif stimPulseAmps(jj) && patternsPlayed(jj) ~= 2*(nOffsets*nPrepulseAmpCombos + 1)
        error(['unexpected number of applications of pattern ' num2str(jj) '... aborting'])
    elseif pSimAmps(jj) && patternsPlayed(jj) ~= 2;
        error(['unexpected number of applications of pattern ' num2str(jj) '... aborting'])
    end
end

%%

%orders movies and times, and checks that all times are within nSamples
for i = 1:length(pattern) %loops through movie chunks
    [times{i} patternOrder] = sort(times{i}, 2, 'ascend'); %#ok<AGROW>
    pattern{i} = pattern{i}(patternOrder); %#ok<AGROW>
    
    if times{i}(end) + 100 > nSamples
        error('at least one of the movie chunks is too long to fit within nSamples')
    end
end

% makes a movie chunk for each row of Pattern and concatenates MovieChunks
MovieChunks = nChunks;
for i = 1:nChunks
    patternChunk = pattern{i};
    timesChunk = times{i};
    
    Chunk = NS_MovieChunkGenerationForExperiment(timesChunk, nSamples, patternChunk);
    MovieChunks = [MovieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end
    
MovieChunksFile=MovieChunks;