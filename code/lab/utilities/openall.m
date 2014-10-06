function openall(mfilename)
% OPENALL     Open all versions of m-file that exist on path (including shadowed versions)
%
% usage:  openall load_params
%
% arguments:  mfilename - The name of the m-file(s) to open
%
% 2010-01 phli

error(nargchk(1, 1, nargin, 'struct'));

% Get all versions of m-file on path, including shadowed versions
mfilepaths = which(mfilename, '-all');

for i = 1:numel(mfilepaths)
    disp(mfilepaths{i});
    open(mfilepaths{i});
end