function average_movie = compute_average_movie_old(spikes,image_spec,varargin)
% compute_average_movie     compute average movie accompanying datarun
%
% usage:  average_movie = compute_average_movie(spikes,image_spec,<params>)
%
% arguments:     spikes - vector of spike times (sec)
%            image_spec - struct with fields:
%
%           	         image_spec.path                path to folder with images
%                        image_spec.data(N).name        name of image stack (relative to the path)
%                                          .index       which image in the stack
%                                          .time        time of image acquisition (sec)
%
%             <params> - struct or list of optional parameters (see below)
%
%
% outputs:  average_movie - struct with fields:
%
%          average_movie.frames      - YxXxT matrix
%                       .frame_count - T-length vector, # images per frame
%                       .avg_frame   - YxX matrix, average of all images (unweighted by spikes)
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
%




% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('subsample', 100);
p.addParamValue('frame_offset', 0.01);
p.addParamValue('range', [-50 50]);
p.addParamValue('normalize','di/i', @(x)any(strcmpi(x,{'di/i','none'})));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% expand params
frame_offset = params.frame_offset;
range = params.range;



% note number of images
num_images = length(image_spec.data);


% compute average

tic


% load one image
temp = load_one_image(image_spec,1);

% initialize movie storage
avg_movie = zeros(size(temp,1),size(temp,2),diff(range));
avg_frame = zeros(size(temp));
frame_count = zeros(diff(range),1);

% initialize storage of all images
all_ims = cell(num_images,1);

% normalization parameter
norm_radius = 40;

% initialize subsampled movie
sub_movie = zeros(size(temp,1),size(temp,2),ceil(num_images/params.subsample));
sub_frame = 1;

% go through each image
for ii = 1:num_images
    
    
    % get the image and time
    [im,im_time] = load_one_image(image_spec,ii);
    
    

    % add it to the overall average frame
    avg_frame = avg_frame + im;
    
    
    % add to subsampled movie
    if mod(ii,params.subsample) == 1
        sub_movie(:,:,sub_frame) = im;
        sub_frame = sub_frame + 1;
    end
    

    % normalize it
    switch params.normalize

        case 'di/i'
            % identify nearby images
            norm_ims = [ii - norm_radius  ii + norm_radius];
            if norm_ims(1) < 1; norm_ims = norm_ims - norm_ims(1) + 1; end
            if norm_ims(2) > num_images; norm_ims = norm_ims - (norm_ims(2)-num_images); end

            % load the sum of nearby images
            norm_sum = zeros(size(temp));
            for jj = norm_ims(1):norm_ims(2)
                % get image jj
                if isempty(all_ims{jj})
                    all_ims{jj} = load_one_image(image_spec,jj);
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
    delta_ts = im_time - spikes;

    % identify which spikes it's close to
    close_spikes = find (   (delta_ts > ( (frame_offset) * range(1) ))  &  (delta_ts < ( (frame_offset) * range(2) ))  );
    
    
    % show progress
    if mod(ii,100)==0
        fprintf('loaded %d of %d in %0.2f sec\n',ii,num_images,toc)
    end
    
    
    % if no spikes are near this image, ignore it
    if isempty(close_spikes)
        continue
    end

    % for each close spike...
    for ss = close_spikes'

        % get offset
        offset = delta_ts(ss);

        % figure out which frame to add it to
        frame = floor(offset/frame_offset) - range(1) + 1;

        %fprintf('spike time: %0.5f    image time: %0.5f    frame number: %d\n',spikes(ss)*1000,im_time*1000,frame)

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

% compute average frame
avg_frame = avg_frame/num_images;



% package output
average_movie.frames = avg_movie;
average_movie.frame_count = frame_count;
average_movie.avg_frame = avg_frame;
average_movie.frame_offset = frame_offset;
average_movie.range = range;
average_movie.sub_movie = sub_movie;



function [im,im_time] = load_one_image(image_spec,index)
% get image and its time

im = double(imread([image_spec.path image_spec.data(index).name],image_spec.data(index).index));

im_time = image_spec.data(index).time;



