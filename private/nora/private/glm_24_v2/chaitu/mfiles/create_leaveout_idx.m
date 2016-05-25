function idx = create_leaveout_idx(tstim,chunksize,sepsize,nchunks,maxt)

% leaving out data from training (originially for cross validation; can be used to avoid training on the boundaries in the block-design experiment)

% tstim:     time duration of single frame [sec]
% chunksize: leaveout chunk size [sec]


% example values in fitting_script.m
% chunksize = 1; [sec]
% sepsize  = 1;
% nchunks = 30;

nframes = floor(chunksize/tstim);
sepframes = floor(sepsize/tstim);

if (sepsize > 0)
    idx = repmat((1:nframes)',1,nchunks) + repmat(makeaxis(sepframes,nchunks),nframes,1);
else
    idx = repmat((1:nframes)',1,nchunks);
end
    
idx = idx(:);
idx = idx(idx < maxt);
