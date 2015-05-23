datarun = load_data('2014-11-24-3/data009');
datarun = load_neurons(datarun);

cell = 5;
spikes = datarun.spikes{cell};
STA = zeros(320,160,30);

% find block start times
block_starts = datarun.triggers([true; diff(datarun.triggers) > 0.9]);

% Load Movie
% movie = get_rawmovie(moviefile, frames) for rawMovie
fitmovie = zeros(320, 160, 120*60);
chunk = 30; % chunks 1 - 30 used for raster
for i_block = 1:10 % the fitting blocks
    
    % Load the movie for that part
    for i_sec = 1:60 % length of fitting segment in seconds
        chunk = chunk+1;
        load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/movie_chunk_' num2str(chunk) '.mat'])
        fitmovie(:,:,120*(i_sec-1)+(1:120)) = movie; 
    end
    
    % Get spike times within the block
    block_start = block_starts(i_block*2);
    block_end = block_starts(i_block*2+1);
    spike_in_block = spikes(spikes > block_start & spikes < block_end) - block_start;
    frame = ceil(spike_in_block*120);
    
    % Calculate STA
    for i_sp = 1:length(frame)
        if frame(i_sp) > 29 && frame(i_sp) < 120*60; STA = STA + fitmovie(:,:,(frame(i_sp)-29):frame(i_sp)); end
    end
    clear spike_in_block

end

% Take a look
for i = 1:30
    imagesc(STA(:,:,i))
    colormap gray
    axis image
    pause(0.1)
end
    