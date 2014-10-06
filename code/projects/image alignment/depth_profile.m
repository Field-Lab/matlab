function dp = depth_profile(I,roi, varargin)
% depth_profile     plot color intensity vs depth in a given roi
%
% usage:  dp = depth_profile(I,F, <params>)
%
% arguments:     I - YxXxCxZ matrix
%              roi - YxX matrix specifying roi, values must be between 0 and 1
%         varargin - struct or list of optional parameters (see below)
%
% outputs:      dp - ZxC matrix of color intensity at each depth
%
%
% optional params, their default values, and what they specify:
%
% norm              'none'       	'how to normalize the curves in each channel
%                                       'none' - don't
%                                       'max' - normalize to have max of 1
%
%
% 2009-09  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('norm','none', @(x)any(strcmpi(x,{'none','max','max,min'})));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% pare down to region of interest

% identify bounds of ROI > 0
[ii jj] = find(roi>0);
roi_x = min(jj):max(jj);
roi_y = min(ii):max(ii);

% pare
I = I(roi_y,roi_x,:,:);
roi = roi(roi_y,roi_x);


% ensure all are double
I = double(I);
roi = double(full(roi));



switch 1

    case 1 % matrixified
        
        % note size
        Is = size(I);

        % reshape image so first dimension is color-depth, second is Y-X
        I = reshape(permute(I,[3 4 1 2]),Is(3)*Is(4),Is(1)*Is(2));

        % reshape roi
        roi = reshape(roi,[],1);

        % multiply
        dp = I*roi;

        % reshape
        dp = reshape(dp,Is(3),Is(4))';

        
    case 2 % loop
        
        dp = zeros(size(I,4),size(I,3));

        for zz = 1:size(I,4)
            for cc = 1:size(I,3)
                dp(zz,cc) = sum(sum(I(:,:,cc,zz).*roi));
            end
        end
        
end


switch params.norm
    case 'max'
        % divide by max
        dp = dp ./ repmat(max(dp,[],1),[size(dp,1) 1]);
    case 'max,min'
        % subtract min
        dp = dp - repmat(min(dp,[],1),[size(dp,1) 1]);
        % divide by max
        dp = dp ./ repmat(max(dp,[],1),[size(dp,1) 1]);
    case 'none'
    otherwise
        error('norm type ''%s'' not recognized',params.norm)
end

