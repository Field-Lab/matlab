function summary_scaled = matrix_scaled_up(summary, scale, params)
% matrix_scaled_up     expand size of matrix by duplicating entries locally
%
%   e.g.  this matrix
%
%       [ 1   5
%         8   2 ]
%
%   becomes the following for a scale factor of 2:
%
%       [ 1   1   5   5
%         1   1   5   5
%         8   8   2   2
%         8   8   2   2 ]
%
%
%  for a 3D matrix, the operation is performed separately for each color
%
%
% usage:  summary_scaled = matrix_scaled_up(summary, scale, params)
%
% arguments:  summary - 3D image matrix (third dimension is color)
%               scale - scale factor (assumed to be the same in x and y, but see below)
%              params - struct of optional parameters (see below)
%
% outputs:
%      summary_scaled - scaled up matrix
%
%
% optional fields in params, their default values, and what they specify:
%
% scale_x, scale_y      <scale>     allows specifying x and y scale independently
%
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.scale_x = scale;
defaults.scale_y = scale;

% combine user and default parameters
params = default_params(defaults, params);


% BODY OF THE FUNCTION

n_x = params.scale_x;
n_y = params.scale_y;

for cc = 1:size(summary,3)
    
    % get this color
    a = summary(:,:,cc);
    
    % scale in x
    b=reshape(repmat(reshape(a,[],1),1,n_y)',n_y*size(a,1),size(a,2));
    
    % scale in y
    c=reshape(repmat(reshape(b',[],1),1,n_x)',n_x*size(b,2),size(b,1))';
    
    % place in output
    summary_scaled(:,:,cc) = c;
end

