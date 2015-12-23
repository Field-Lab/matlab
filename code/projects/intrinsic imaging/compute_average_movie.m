function average_movie = compute_average_movie(triggers, stack, varargin)
% compute_average_movie     compute average movie accompanying datarun
%
% usage:  average_movie = compute_average_movie(triggers,stack,<params>)
%
% arguments: triggers - vector of spike times (sec)
%            stack - struct with fields:
%
%            stack.path                path to folder with images
%            stack.data(N).name        name of image stack (relative to the path)
%                         .index       which image in the stack
%                         .time        time of image acquisition (sec)
%
%            <params> - struct or list of optional parameters (see below)
%
%
% outputs: average_movie - struct with fields:
%
%          average_movie.frames      - YxXxT matrix
%                       .frame_count - T-length vector, # images per frame
%                       .average_frame - Average of the triggered frames
%                       .full_average   - YxX matrix, average of all images (unweighted by triggers)
%                       .sub_movie   - subset of the original movie
%                       .frame_offset - value of input parameter
%                       .range       - value of input parameter
%
%
%
%
% optional params, their default values, and what they specify:
%
% verbose           false         	show output
% subsample         []              period of frames put into the subsampled movie
%                                       if empty, don't make subsampled movie
% frame_offset      0.01            size of time bin in which to combine images (sec)
% range             [-50 50]        range of time bins around the spike
%
%
% 2010-03  gauthier
% 2010-07  phli - ported to new stack
%


% SET UP OPTIONAL ARGUMENTS
p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('subsample', 100);
p.addParamValue('frame_offset', 0.01);
p.addParamValue('range', [-50 50]);
p.addParamValue('normalize','none', @(x)any(strcmpi(x,{'di/i','none'})));
p.addParamValue('norm_radius', 40);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% expand params
frame_offset = params.frame_offset;
range = params.range;

% note number of images
num_images = stack_length(stack);

% compute average
tic

% load info for images, use first image info as template for malloc
stack = load_slice_imfs(stack);
imf = slice_imf(stack, 1);

% initialize movie storage
avg_movie = zeros(imf.Height, imf.Width, diff(range));
full_average = zeros(imf.Height, imf.Width);
frame_count = zeros(diff(range), 1);

% initialize storage of all images
all_ims = cell(num_images, 1);

% initialize subsampled movie
sub_movie = zeros(imf.Height, imf.Width, ceil(num_images/params.subsample));
sub_frame = 1;

% go through each image
for ii = 1:num_images
    
    % get the image data, info, and time
    im = get_slice(stack, ii);
    im = double(im);
    imf = slice_imf(stack, ii);
    im_time = get_slice_trigger(stack, ii);
    

    % add it to the overall average frame
    full_average = full_average + im;
    
    
    % add to subsampled movie
    if mod(ii,params.subsample) == 1
        sub_movie(:,:,sub_frame) = im;
        sub_frame = sub_frame + 1;
    end
    

    % normalize it
    switch params.normalize

        case 'di/i'
            % identify nearby images
            norm_radius = params.norm_radius;
            norm_ims = [ii - norm_radius  ii + norm_radius];
            if norm_ims(1) < 1; norm_ims = norm_ims - norm_ims(1) + 1; end
            if norm_ims(2) > num_images; norm_ims = norm_ims - (norm_ims(2)-num_images); end

            % load the sum of nearby images
            norm_sum = zeros(imf.Height, imf.Width);
            for jj = norm_ims(1):norm_ims(2)
                % get image jj
                if isempty(all_ims{jj})
                    all_ims{jj} = double(get_slice(stack,jj));
                    %fprintf('loaded %d\n',jj)
                end
                % add it to sum
                norm_sum = norm_sum + all_ims{jj};
            end
            norm_avg = norm_sum / (diff(norm_ims) + 1);

            % delete old images
            if norm_ims(1) > 1
                all_ims{norm_ims(1)-1} = [];
                %fprintf('   deleted %d\n',norm_ims(1)-1)
            end

            %fprintf('%d: %d bytes\n',ii,getfield(whos('all_ims'),'bytes'))

            % use to normalize image
            im = (im - norm_avg) ./ norm_avg;
            
        case 'none'
            
        otherwise
            error('normalization type ''%s'' not recognized',params.normalize)
    end



    % compute offset from each spike
    delta_ts = im_time - triggers;

    % identify which triggers it's close to
    close_triggers = find (   (delta_ts > ( (frame_offset) * range(1) ))  &  (delta_ts < ( (frame_offset) * range(2) ))  );
    
    
    % show progress
    if params.verbose && mod(ii,100) == 0
        fprintf('loaded %d of %d in %0.2f sec\n', ii, num_images, toc)
    end

    % for each close spike...
    for ss = close_triggers'

        % get offset
        offset = delta_ts(ss);

        % figure out which frame to add it to
        frame = floor(offset/frame_offset) - range(1) + 1;

        %fprintf('spike time: %0.5f    image time: %0.5f    frame number: %d\n',triggers(ss)*1000,im_time*1000,frame)

        % add it
        avg_movie(:,:,frame) = avg_movie(:,:,frame) + im;
        frame_count(frame) = frame_count(frame) + 1;

    end

end

% divide each frame by number of images
for ff=1:length(frame_count)
    if frame_count(ff) ~= 0
        avg_movie(:,:,ff) = avg_movie(:,:,ff) / frame_count(ff);
    end
end

% compute full average
full_average = full_average/num_images;

% Compute average of triggered frames
average_frame = mean(avg_movie, 3);

% package output
average_movie.frames = avg_movie;
average_movie.frame_count = frame_count;
average_movie.full_average = full_average;
average_movie.average_frame = average_frame;
average_movie.frame_offset = frame_offset;
average_movie.range = range;
average_movie.sub_movie = sub_movie;