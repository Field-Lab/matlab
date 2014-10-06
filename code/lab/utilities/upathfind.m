function [found, matches] = upathfind(glob, pathselect)
% UPATHFIND     Glob search user path for matching files
% usage: [found, matches] = upathfind(glob, pathselect)
%
% Searches the Matlab path for files matching pattern GLOB.  Excludes
% Matlab's internal functions and toolbox functions.  Something similar to 
%   !find MATLAB-STANDARD -name GLOB
%
% 2012-07 phli
%

% Get all path entries as cell array of strings
p = split(path, ':')';

% Filter out paths matching matlabroot, as we are looking for things only on the user path
p = select(p, @userpath_select);

% Process additional path selection?
if nargin > 1 
    if ischar(pathselect)
        % Convert pathselect string into selection function
        pathselect = @(str) (~isempty(strfind(str, pathselect)));
    end
    p = select(p, pathselect);
end

% Search remaining paths for matching strings
matches = cellfun(@(d)(dir(sprintf('%s/%s', d, glob))), p, 'UniformOutput', false);
matches = cell2mat(matches);
found = {matches.name}';

if nargout == 0
    cellfun(@(str)(disp(['    ' str])), found);
    clear found;
end


function is_userpath = userpath_select(path)
is_userpath = isempty(strfind(path, matlabroot()));