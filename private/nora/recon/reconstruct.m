function movie = reconstruct(spikes, linear_filter, varargin)

% linear filter should be in [space, space, time]
% spikes should be in seconds

p = inputParser;
p.addParameter('fps', 1/0.0083275, @isnumeric) % frames per second, using calculated tstim 
p.parse(varargin{:});

% round to nearest whole second
movie_length = ceil(ceil(max(spikes)) * p.Results.fps);

% these are fixed for this function
filter_size = size(linear_filter);
movie = zeros([filter_size(1:2), movie_length]);

% weird indexing is to take care of spikes in the beginning frames
for i_spike = 1:length(spikes)
   spike_frame = floor(spikes(i_spike) * p.Results.fps);
   idx_begin = (spike_frame-filter_size(3)+1);
   if idx_begin < 1
       idx_begin = 1; 
   end
   movie(:,:,idx_begin:spike_frame) = movie(:,:,idx_begin:spike_frame) + linear_filter(:,:,(end-(spike_frame-idx_begin)):end);
end

end