function path2data=find_data(piece, run, streamed, online)

if online % during acquisition
    path2data=fullfile('/Volumes/Acquisition/Analysis/', piece, run, run);
elseif streamed % offline, streamed data
    path2data=fullfile('/Volumes/Analysis/', piece, 'streamed', run, run);
else % offline, full data
    path2data=fullfile('/Volumes/Analysis/', piece, run, run);
end