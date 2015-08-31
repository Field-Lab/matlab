function movie = lpf_test(varargin)
%% Parse optional input 
p = inputParser;
p.addParamValue('tau', 2)
p.addParamValue('sigmas', 6)
p.addParameter('sigma', 2)
p.addParameter('dim', [320 640])
p.addParameter('seed', 11111)
p.addParameter('rng-type', 'java')
p.addParameter('back-rgb', 0.5)
p.addParameter('frames', 900)
p.parse(varargin{:});

%% Make the filters

% Scalars for adding time
scalar_previous = exp(-1/p.tau);
scalar_current = sqrt(scalar.previous); % not sure what this is supposed to be

% Gaussian spatial filter
gauss_filter = make_gaussian_filter(p.sigma, p.sigmas);

% Gaussian lookup table
% magic numbers haha
lookup_table = icdf('Normal', (1:1024)/2025, 0.5, 0.16);

%% Make the movie!

% Initialize
frame_previous = zeros(p.dim);
movie = zeros([p.frames, dim]);

% Make frame by frame
for i_frame = 1:p.frames
    
    % keeping tabs on the progress
    if ~mod(i_frame, 100)
        disp(['Frame ' num2str(i_frame)])
    end
    
    % make frame
    frame_current = randn(p.dim); % replace with java rand
    frame_current = lookup_table(frame_current);
    
    % filter
    movie(i_frame, :, :) = conv2(scalar_current*frame_current + scalar_previous*frame_previous, gauss_filter);

    % save for next frame
    frame_previous = frame_current;
   
end

% round to 0 to 255
movie = movie * 256;
movie = floor(movie);
if min(movie(:)) < 0 || max(movie(:)) > 255
    error('Movie out of range')
end

end

function filter = make_gaussian_filter(sigma, sigmas)

% calculate the filter size in number of pixels
filter_size = round(sigma*sigmas);

% make it odd
if ~mod(filter_size, 2)
    filter_size = filter_size+1;
end

% Make the 1D window
oneD_filter = gausswin(filter_size, sigma);

% Make the 2D filter
filter = oneD_filter'*oneD_filter;

% Normalize the filter
filter = filter/norm(filter(:)); %potentially need to add 0.16 factor?

end


