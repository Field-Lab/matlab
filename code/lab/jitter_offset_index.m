function offsets = jitter_offset_index(n_x,n_y)
% jitter_offset_index     create offset for jitter stimulus analysis
%
% usage:  offsets = jitter_offset_index(n_x,n_y)
%
% arguments:      n_x - number of jitter steps in x
%                 n_y - in y
%
% 	if only one argument is specified, it is used for n_x and n_y
%
% outputs:    offsets - result of computation
%
%
% 2008-12 gauthier
%

% ensure n_y exists
if ~exist('n_y','var')
    n_y = n_x;
end

% x offsets
temp_x = repmat(0:n_x-1,n_x,1);
offsets_x = reshape(temp_x,1,[]);

% y offsets
temp_y = repmat(0:n_y-1,n_y,1);
offsets_y = reshape(temp_y',1,[]);

% combine
offsets = [offsets_x - floor(n_x/2); offsets_y - floor(n_y/2)];

