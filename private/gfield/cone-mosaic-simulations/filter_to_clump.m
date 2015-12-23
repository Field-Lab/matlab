function [cone_types, extras] = filter_to_clump(locations, filter_sigma, pL, delta_pL, varargin)
%
% spatially filter a cone mosaic to produce a tendency of L an M cones to
% clump.
%
% usage:     cone_types = filter_to_clumping(dim_size, buffer, scale, pL, deltapL)
%
%  The seed for this code was provided by Hie


% PARSE INPUT
p = inputParser;
p.addRequired('locations', @isnumeric)
p.addRequired('filter_sigma', @isnumeric)
p.addRequired('pL', @isnumeric)
p.addRequired('delta_pL', @isnumeric)

p.addParamValue('buffer', 15, @isnumeric) % spatial buffer for removing edge effects of filtering
p.addParamValue('seed', [], @isnumeric)
p.addParamValue('filter_size', 4, @isnumeric) % radius of filter in units of sigma

p.parse(locations, filter_sigma, pL, delta_pL, varargin{:})

% INITIALIZE VARIABLES
buffer = p.Results.buffer;
seed = p.Results.seed;
filter_size = p.Results.filter_size;

% seed random number generated
if ~isempty(seed)
    rand('twister', seed)
end

% BODY OF FUNCTION
base_size = ceil(max(max(locations))); % the +1 insures the profile field is large enough to capture all cones
%base_size = sqrt(length(locations(:,1))); % assuming a square area of points this provides the dimensions
dim_size = base_size + 2 * buffer; % adds buffer (padding) which will be removed after filtering

% generate a white noise from uniform distribution
profile = rand(dim_size);

%now filter so that there will be clumps with spatial scale 'scale'
f = fspecial('gaussian', [filter_size*filter_sigma], filter_sigma);
profile = conv2(profile,f,'same');

% remove buffer
profile = profile(buffer+1:(dim_size-buffer), buffer+1:(dim_size-buffer));

%now rescale the values in profile to fill the desired range along with
% ad hoc smoothing
av = mean(mean(profile));
maxp = mean(max(profile));
minp = mean(min(profile));
for i=1:base_size
    for j=1:base_size
        if profile(i,j)>= maxp
            profile(i,j)=maxp;
        else if profile(i,j) <=minp
                profile(i,j)=minp;
            end
        end
    end
end
range = maxp - minp;
profile = profile - (av .* ones(base_size));

% make profile mean about 0.5 and max and min about 0 and 1
profile = av + (profile ./ range); 

%rescale for close to the desired pL and delta pL (actual values for each
%mosaic will be returned
if delta_pL < 1-(2*abs(0.5-pL))
    profile = profile - 0.5 * ones(base_size,base_size);
    profile = pL + profile .* delta_pL;
else if delta_pL >= 1-(2*abs(0.5-pL))
    profile = profile - 0.5*ones(base_size,base_size);
    profile = pL+profile .* 2 * (delta_pL - 0.5 + abs(0.5 - pL));  
    end
end

% ensure probabilities are withing limits of 0 and 1
for i=1:base_size
    for j=1:base_size
        if profile(i,j)>1
            profile(i,j)=1;
        else
            if profile (i,j)<0;
                profile(i,j)=0;
            end
        end
    end
end


% now assign cones based on these local probabilities 0=L 1=M
for i=1:length(locations(:,1))
    temp_indices = ceil(locations(i,:)); 
    if rand(1) > 1 - profile(temp_indices(1), temp_indices(2))
      cone_types(i) = 'L';
      used_pL(i) = profile(temp_indices(1), temp_indices(2));
    else
      cone_types(i) = 'M';
      used_pL(i) = profile(temp_indices(1), temp_indices(2));
    end
end

l_cone_ids = find(cone_types == 'L');
m_cone_ids = find(cone_types == 'M');

%now get actual average pL and range pL values
actual_pL = mean(used_pL);
actual_delta_pL = max(max(used_pL)) - min(min(used_pL));

% pack additional information into 'extras' struct for return
extras.l_cone_indices = l_cone_ids;
extras.m_cone_indices = m_cone_ids;
extras.actual_pL = actual_pL;
extras.actual_delta_pL = actual_delta_pL;
