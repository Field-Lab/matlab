function out_mat = matrix_jitter_subsample(matrix,n,offset,varargin)
% jitter_subsample     subsamples and jitters a matrix
%
%   see documentation file jitter/jitter_subsample for details
%
%
% usage:  out_mat = matrix_jitter_subsample(matrix,n,offset,varargin)
%
% arguments:   matrix - 3D matrix
%                   n - pixels per stixel
%              offset - length 2 vector, first entry is x offset, second is y
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    out_mat - result
%
%
% optional params, their default values, and what they specify:
%
% aperture_x     	':'     which x values to return.  if ':', return all
% aperture_y        ':'     same for y
% mask_value        0.5     values to return in the mask region
%
%
% 2008-12 gauthier
%


% SET UP OPTIONAL ARGUMENTS

% specity required parameters
q = inputParser;
q.addRequired('matrix',@(x) min(size(x,1),size(x,2))>=3  );

% verify them
q.parse(matrix);


% specify list of optional parameters
p = inputParser;
p.addParamValue('aperture_x', ':', @(x)length(x)==length(intersect(1:size(matrix,2)*n,x)));
p.addParamValue('aperture_y', ':', @(x)length(x)==length(intersect(1:size(matrix,1)*n,x)));
p.addParamValue('mask_value', 0.5,@(x)numel(x)==1);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION

% scale up matrix
big_mat = matrix_scaled_up(matrix,n);

% initialize output matrix
out_x = size(big_mat,2);
out_y = size(big_mat,1);
out_mat = params.mask_value*ones(out_y,out_x,size(matrix,3));

% insert values of scaled up original matrix
out_mat(n+1:end-n,n+1:end-n,:) = big_mat(n+1-offset(2):end-n-offset(2),n+1-offset(1):end-n-offset(1),:);

% only return aperture, is specified
if ~strcmpi(params.aperture_x,':')
    out_mat = out_mat(:,params.aperture_x,:);
end
if ~strcmpi(params.aperture_y,':')
    out_mat = out_mat(params.aperture_y,:,:);
end



