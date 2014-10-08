function MovieChunksFile = generateMovieSpatialPatternsVary1Elec(nElec, periodInMs, offsets, orders, variedElec)

% generates movieChunksFile for stimulus file generation
%
% offsets: a vector of offset values specifying delays between subsequent pulses within a sequence
% (in ms)
%
% specifies order in which patterns are played during stimulus (time each pattern is played)
%
% movie chunk length is limited to <=2 seconds, so patterns must be broken up into multiple movie chunks, corresponding
% to different rows in Pattern


nOffsets = length(offsets);
nOrders = size(orders, 1);

period = periodInMs*20;
offsets = offsets*20;


% length of each movie chunk (samples)
nStimuliPerElecAmp = nOrders*nOffsets + 2; %sequences + simultaneous + varied electrode alone


%figure out how to break up into movie chunks
maxPatternsPerChunk = floor(40000/period);
nChunksPerElecAmp = ceil(nStimuliPerElecAmp/maxPatternsPerChunk);
stimPerChunk = ceil(nStimuliPerElecAmp/nChunksPerElecAmp);

nChunksTotal = nChunksPerElecAmp*nElec*5; %5 amplitudes to test for each electrode in sequence

nSamples = ceil((stimPerChunk*period)/2000)*2000; %round to nearest tenth of a second
h = msgbox([num2str(nChunksTotal) ' movie chunks required', 10, 'Set trigger interval for ''spatial_vary'' to ' num2str(nSamples/20000) ' seconds']);
boxPos = get(h, 'position');
set(h, 'position', [400 500 boxPos(3) boxPos(4)])




%%
  
TimesAll = cell(1, nChunksTotal);
PatternAll = cell(1, nChunksTotal);

iChunk = 0;
for ii = 1:nElec %electrode that's varying
    varyPatterns = find(variedElec == ii);
    for jj = 1:5 %which amplitude
        
        %simultaneous
        Pattern = varyPatterns(jj);
        Times = 0;
        
        patternsInSequences = 1:nElec; %first nElec patterns are full amplitudes
        patternsInSequences(ii) = varyPatterns(jj+5); %replace varied electrode with near-thresh amp
        
        %sequences
        for mm = 1:nOffsets
            for kk = 1:nOrders
                patternTime = period*((mm-1)*nOrders + kk); %start time for this offset/order
                
                Pattern = [Pattern patternsInSequences(orders(kk,:))];
                
                tempTimes = 0:offsets(mm):offsets(mm)*(nElec-1);
                Times = [Times tempTimes + patternTime];
            end
        end
        
        %electrode alone
        Pattern = [Pattern varyPatterns(jj+5)];
        Times = [Times period*(nOrders*nOffsets+1)];

        
        %break into chunks
        for kk = 1:nChunksPerElecAmp
            iChunk = iChunk + 1;
            stimThisChunk = (Times < stimPerChunk*period);
            TimesAll{iChunk} = Times(stimThisChunk);
            PatternAll{iChunk} = Pattern(stimThisChunk);

            %what's left
            Times = Times(~stimThisChunk);
            if kk < nChunksPerElecAmp
                Times = Times - Times(1); %restart at time = 0
            end
            Pattern = Pattern(~stimThisChunk);
        end
        if ~isempty(Times) || ~isempty(Pattern)
            errordlg('something went wrong when breaking movie into chunks')
        end
    end
end


% generate movie chunks
MovieChunksFile = nChunksTotal;

%patternsUsed = zeros(1, length(variedElec));
for ii = 1:nChunksTotal
%     %sanity check
%     for jj = 1:length(patternsUsed)
%         patternsUsed(jj) = patternsUsed(jj) + sum(PatternAll{ii} == jj);
%     end
    
    Chunk = NS_MovieChunkGenerationForExperiment(TimesAll{ii}, nSamples, PatternAll{ii});
    MovieChunksFile = [MovieChunksFile Chunk];
end

