function [movie, p] = lpf_make_mat(matfile_save_location, varargin)
% NB Assuming BW for now 2015-08-17

%% Parse optional input 
p = inputParser;
% p.addRequired('matfile_save_location'); % where to save the matfiles
p.addParameter('tau', 12) % sets scale of decay over time
p.addParameter('sigmas', 6) % sets size of filter relative to gaussian width
p.addParameter('sigma', 4) % sets size of gaussian
p.addParameter('dim', [320 160]) % size of screen 
p.addParameter('rng_type', 'java') % random number generator to use
p.addParameter('back_rgb', 0.5)  % back rgb
p.addParameter('seconds', 1) % number of seconds at 120Hz to generate
p.addParameter('independent', 0); % BW or RGB
p.addParameter('intensity_SD', 0.16); % SD of pixels
p.addParameter('seed', 11111);
p.parse(varargin{:});

%% Make the filters
 
% Scalars for adding time
scalar_previous = exp(-1/p.Results.tau);
scalar_current = sqrt(1 - scalar_previous^2);

% Gaussian spatial filter

% calculate the filter size in number of pixels
filter_size = round(p.Results.sigma*p.Results.sigmas);

% make it odd
if ~mod(filter_size, 2)
    filter_size = filter_size+1;
end

% Make the 1D window
alpha = (filter_size-1)/(2*p.Results.sigma); % just from matlab's weird code definition for gausswin
oneD_filter = gausswin(filter_size, alpha); % don't need the two pi because just normalize later 

% Make the 2D filter
gauss_filter = oneD_filter*oneD_filter';

% Normalize the filter
gauss_filter = gauss_filter/norm(gauss_filter(:)); % removed 0.16 factor here because it is in the frame function from alex

clear oneD_filter alpha filter_size

%% Make the movie!

% Initialize

% frame_previous = zeros(p.Results.dim);
% movie = zeros([p.Results.dim, p.Results.frames], 'uint8'); % check where frames go? also should this be times 3? I think so 
movie = zeros([p.Results.dim, 120]); % make a one second movie
state = Init_RNG_JavaStyle(p.Results.seed, p.Results.dim, p.Results.back_rgb);

% lookup table size
n_table = 1024;

% get stuff ready to pass to movie calling function. All this is stuff from
% Alex... not entirely sure what it means
rgb_vect = [1 1 1]*p.Results.intensity_SD; % for gaussian
back_rgb = round([1 1 1]*p.Results.back_rgb*255);
if p.Results.independent % RGB
    noise_type = 3;
else % BW
    noise_type = 2;
end
n_bits = 16; %8;
tmp = repmat(norminv((1:n_table)/(n_table+1), 0, 1)', 65536/n_table, 1);
tmp = repmat(tmp,1,3) .* repmat(rgb_vect,65536,1);%  + repmat(back_rgb, 65536, 1);
% tmp = uint8(round(255 * tmp))';
tmp = tmp';
lut = single(tmp(:));

% Make one second at a time
for i_sec = 1:p.Results.seconds
    
    disp(i_sec)
    
    % Make frames for that second
    for i_frame = 1:120
        
        frame = Draw_Random_Frame_float(state, p.Results.dim(1), p.Results.dim(2), lut, ...
            [], uint8(back_rgb), p.Results.dim(1), p.Results.dim(2), noise_type, n_bits, 1);
        frame_current = double(squeeze(frame(1,:,:)));
        
        % filter in time
        if i_frame > 1 || i_sec > 1
            unfiltered_frame = scalar_current*frame_current + scalar_previous*frame_previous;
        else
            unfiltered_frame = frame_current;
        end
        frame_previous = unfiltered_frame;
        
        % filter in space
        filtered_frame = imfilter(unfiltered_frame, gauss_filter, 'circular');
        
        % save to movie array
        movie(:, :, i_frame) = floor(256*(filtered_frame+p.Results.back_rgb));
        
    end
    
    % round to 0 to 255
    movie(movie > 255) = 255;
    movie(movie < 0) = 0;
    if min(movie(:)) < 0 || max(movie(:)) > 255
        error('Movie out of range')
    end
    
    % Save the movie chunk to matfile
    save([matfile_save_location 'movie_chunk_' num2str(i_sec) '.mat'], 'movie')
    
end

% Save the movie parameters to file
save([matfile_save_location 'movie_params.mat'], 'p')




end