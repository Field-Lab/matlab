function [ta, e_ta] = compute_ta_from_ind(tai, imfolder, zeroval)
%COMPUTE_TA_FROM_IND Computes a triggered average from a matrix of indices
%generated by compute_ta_ind and an events folder.
%
%  TA = COMPUTE_TA_FROM_IND(TAI, IMFOLDER) will use the events stored in
%  the folder IMFOLDER (assumed to be stored by chunks in alphanumerical
%  older, and saved as mat files) and the matrix of indices TAI to
%  calculate a spike-triggered average of the events.
%  
%  TA = COMPUTE_TA_FROM_IND(..., ZEROVAL) speficies the neutral value for
%  this TA computation. For example, it should be 0.5 for white noise
%  movies, 128 for natural scenes, etc. Defaults to 0.5.
%
%  [~, E_TA] = COMPUTE_TA_FROM_IND further returns the standard error 
%  of the mean for the triggered average.

% Maximum number of movie chunks in memory at one time
MAX_CHUNKS_LOADED = 3; 
movie_chunks = cell(MAX_CHUNKS_LOADED, 1);
chunksloaded = zeros(MAX_CHUNKS_LOADED, 1);

if nargin == 2
    zeroval = 0.5;
end
neventsaveraged = size(tai, 1);
tadepth = size(tai, 2);
nevents = numel(tai);

% We'll work on the transpose of the TA indices matrix, as it is more
% natural to index in column-major mode in Matlab.
tai = tai.';

% First, we need to figure out where all the event chunks are stored, and
% in which chunk we should be looking for each event.

% Get the list of all chunk files
allchunks = dir(fullfile(imfolder, '*.mat'));

% Open the first one to figure out how many events are stored per chunk
load(fullfile(imfolder, allchunks(1).name));
eventsperchunk = length(frames);

% Event 0 is the special event "no frame" (gray), create it here
eventzero = zeroval*ones(size(frames{1}));

% Now we can calculate the TA and its associated error
ta = repmat({zeros(size(frames{1}))}, tadepth, 1);
e_ta = repmat({zeros(size(frames{1}))}, tadepth, 1);
event_counter = zeros(tadepth, 1);

for kk = 1:nevents
    cevent = tai(kk);
    cframe = mod(kk-1, tadepth) + 1;
    
    % Find which chunk corresponds to the event and what its index in the
    % chunk is
    chunkindex = ceil(cevent/eventsperchunk);
    posinchunk = mod(cevent, eventsperchunk);
    if posinchunk == 0
        posinchunk = eventsperchunk;
    end
    
    % Add the event accordingly to mean and variance
    % Event 0 is a special case handled differently
    if tai(kk)>0
        
        % We only load a movie chunk if it wasn't already in memory
        if sum(chunksloaded == chunkindex) == 0
            % Pop the last chunk
            movie_chunks(2:end) = movie_chunks(1:(end-1));
            chunksloaded(2:end) = chunksloaded(1:(end-1));
            
            % Push in the new chunk
            load(fullfile(imfolder, allchunks(chunkindex).name));
            movie_chunks{1} = frames;
            chunksloaded(1) = chunkindex;
            
            % We'll use this chunk
            bufchunkind = 1;
        else
            bufchunkind = find(chunksloaded == chunkindex);
        end
        
        % For each time the event contributed, update mean/variance.
        % This is our new sample
        x = movie_chunks{bufchunkind}{posinchunk} - zeroval;
        % Update event counter
        event_counter(cframe) = event_counter(cframe) + 1;
        % Running mean and variance
        delta = x - ta{cframe};
        ta{cframe} = ta{cframe} + delta/event_counter(cframe);
        e_ta{cframe} = e_ta{cframe} + delta.*(x - ta{cframe});
        
    else
        
        % New sample
        x = eventzero - zeroval;
        % Update event counter
        event_counter(cframe) = event_counter(cframe) + 1;
        % Running mean and variance
        delta = x - ta{cframe};
        ta{cframe} = ta{cframe} + delta/event_counter(cframe);
        e_ta{cframe} = e_ta{cframe} + delta.*(x - ta{cframe});
        
    end
end

% Normalize things
for kk=1:length(e_ta)
    % Go from M2 to variance
    e_ta{kk} = e_ta{kk} / (event_counter(kk) - 1);
    % And we actually want to return the standard error of the mean
    e_ta{kk} = sqrt(e_ta{kk})/sqrt(neventsaveraged);
end

end % compute_ta_from_ind