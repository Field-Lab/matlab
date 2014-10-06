function MovieChunksFile = generateMovieSpatialPatternsVaryAllElecs(nElec, periodInMs, offsets, orders)

% offsets: a vector of offset values specifying delays between subsequent pulses within a sequence
% (in ms)
%
% assumes that first stimulus in Array is simultaneous stimulation and the remaining stimuli are
% each individual electrode
%

nOffsets = length(offsets);
nOrders = size(orders, 1);

period = periodInMs*20;
offsets = offsets*20;


% length of each movie chunk (samples)
nStimuli = nOrders*nOffsets + 1 + nElec; %(sequences + simultaneous + individual electrodes)
%samplesRequired = nStimuli*period;

%figure out how to break up into movie chunks
maxPatternsPerChunk = floor(40000/period);
nChunks = ceil(nStimuli/maxPatternsPerChunk);
stimPerChunk = ceil(nStimuli/nChunks);


nSamples = ceil((stimPerChunk*period)/2000)*2000; %round to nearest tenth of a second
h = msgbox([num2str(nChunks) ' movie chunks required', 10, 'Set trigger interval for ''spatial_short'' to ' num2str(nSamples/20000) ' seconds']);
boxPos = get(h, 'position');
set(h, 'position', [400 500 boxPos(3) boxPos(4)])


%% specify order/timing of pattern application
    
orders = orders + 1; %because the first pattern corresponds to simultaneous stimulation

nPatterns = (nOffsets*nOrders+1)*nElec + 1; %(sequences + individual elecs) + simultaneous

%simultaneous
Pattern = 1;
Times = 0;

%sequences
for ii = 1:nOffsets
    for jj = 1:nOrders
        patternTime = period*((ii-1)*nOrders + jj); %start time for this offset/order

        Pattern = [Pattern orders(jj,:)];

        tempTimes = 0:offsets(ii):offsets(ii)*(nElec-1);
        Times = [Times tempTimes + patternTime];
    end
end

%individual electrodes
Pattern = [Pattern (1:nElec)+1];
Times = [Times period*(nOrders*nOffsets+1):period:period*(nOrders*nOffsets+nElec)];


%% break into movie chunks

MovieChunksFile = nChunks;
for ii = 1:nChunks
    stimThisChunk = (Times < stimPerChunk*period);
    
    Chunk = NS_MovieChunkGenerationForExperiment(Times(stimThisChunk), nSamples, Pattern(stimThisChunk));
    MovieChunksFile = [MovieChunksFile Chunk];
    
    Times = Times(~stimThisChunk);
    if ii < nChunks
        Times = Times - Times(1); %restart at time = 0
    end
    Pattern = Pattern(~stimThisChunk);
end

if ~isempty(Times) || ~isempty(Pattern)
    errordlg('something went wrong when breaking movie into chunks')
end


