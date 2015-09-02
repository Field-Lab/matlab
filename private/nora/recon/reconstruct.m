function movie = reconstruct(spikes, linear_filter, varargin)

% linear filter should be in [space, space, time]
% spikes should be in seconds

% can turn into varargin
frames_per_second = 120;

% round to nearest whole second
movie_length = ceil(max(spikes)) * frames_per_second;

% these are fixed for this function
filter_size = size(linear_filter);
movie = zeros([filter_size(1:2), movie_length]);


% weird indexing is to take care of spikes in the beginning frames
for i_spike = 1:length(spikes)
   spike_frame = floor(spikes(i_spike) * frames_per_second);
   idx_begin = (spike_frame-filter_size(3)+1);
   if idx_begin < 1
       idx_begin = 1; 
   end
   movie(:,:,idx_begin:spike_frame) = movie(:,:,idx_begin:spike_frame) + linear_filter(:,:,(end-(spike_frame-idx_begin)):end);
end


% Not sure if these should be in there
movie = movie - min(movie(:));
movie = movie/max(movie(:));

end