function offsets =  compute_jitter_offsets(num_frames,mdf_seed,stixel_height,stixel_width)
% compute_jitter_offsets     compute what the jitter was in each white noise movie frame
%
% usage:  offsets =  compute_jitter_offsets(num_frames,mdf_seed,stixel_height,stixel_width)
%
% arguments:     num_frames - return offsets for the first num_frames frames
%                                try  ceil(length(datarun.triggers)*100/datarun.stimulus.interval)
%                  mdf_seed - seed of the movie file, usually 11111
%             stixel_height - height of each stixel, in pixels
%              stixel_width - width
%
% outputs:    offsets - Nx2 matrix of jitter offsets.  first column is x, second is y.
%
%
% 2008-12 gauthier
%



% initialize variable to store offsets
rand_nums_x = zeros(num_frames,1);
rand_nums_y = zeros(num_frames,1);

% initialize random number generator
rng = edu.ucsc.neurobiology.vision.math.RandomJavaV2(mdf_seed);

% get seeds
for rr = 1:size(rand_nums_x,1)
    rand_nums_x(rr,1) = rng.nextBits(16);
    rand_nums_y(rr,1) = rng.nextBits(16);
end

% compute offset of each frame
offsets_x = mod(rand_nums_x,stixel_width) - floor(stixel_width/2);
offsets_y = mod(rand_nums_y,stixel_height) - floor(stixel_height/2);

% package into a single matrix
offsets = [offsets_x offsets_y];
