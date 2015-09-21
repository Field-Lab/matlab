function movie = reconstruct(spikes, linear_filter, varargin)

% linear filter should be in cell [space, space, time]
% spikes should be in cells in seconds

p = inputParser;
p.addParameter('fps', 1/0.0083275, @isnumeric) % frames per second, using calculated tstim 
p.parse(varargin{:});

% round to nearest whole second
movie_length = ceil(ceil(max(spikes{2})) * p.Results.fps);

% these are fixed for this function
filter_size = size(linear_filter{1});
movie = zeros([filter_size(1:2), movie_length]);

% cycle through cell spikes
% weird indexing is to take care of spikes in the beginning frames
for i_cell = 1:length(spikes)
    for i_spike = 1:length(spikes{i_cell})
        spike_frame = floor(spikes{i_cell}(i_spike) * p.Results.fps);
        idx_begin = (spike_frame-filter_size(3)+1);
        if idx_begin < 1
            idx_begin = 1;
        end
        movie(:,:,idx_begin:spike_frame) = movie(:,:,idx_begin:spike_frame) + linear_filter{i_cell}(:,:,(end-(spike_frame-idx_begin)):end);
    end
end

end