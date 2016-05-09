for i = 1:100
    load(['/Volumes/Lab/Users/Nora/new_stim_nora/NSEM intervals/matfiles/movie_chunk_' num2str(i) '.mat'])
    interval(i) = size(movie_chunk, 3);
end
default_colors = get(gca,'ColorOrder');
interval_time = cumsum(interval/120);

%%

[~, ~, spikes] = interleaved_data_prep(datarun, [1800, 1800], 20, 'cell_spec', 4461 );
i = 0;
time = 0;
while time < 15
    i = i+1;
    hold on; plot([interval_time(i), interval_time(i)], [0 20], 'Color', default_colors(1,:))
    time = interval_time(i+1);
end


