function gm = fitGMM(binnedResponses,stimulus)

tsp = find(binnedResponses==1);
spike_stim = stimulus(tsp,:);

options = statset('Display','final');
gm = fitgmdist(spike_stim,4,'Options',options);

end