function p = genpath2(d, skips)
%GENPATH2 Generate recursive toolbox path, excluding specified directories
%   P = GENPATH(D, SKIPS) returns a path string starting in D, plus, 
%   recursively, all the subdirectories of D, excluding directories
%   matching SKIPS, a cell array of strings.
%
%   phli 2010-05, copied and modified from original MATLAB genpath
%

% initialise variables
methodsep = '@';  % qualifier for overloaded method directories
p = '';           % path to be returned

% Generate path based on given root directory
files = dir(d);
if isempty(files)
    return
end

% Add d to the path
p = [p d pathsep];

% select only directory entries from the current listing
isdir = logical(cat(1,files.isdir));
dirs = files(isdir);


%
% Recursively descend through directories which are neither
% private nor "class" directories nor SVN directories.
%
for i=1:length(dirs)
    dirname = dirs(i).name;
    
    if any(skips, @(skip)(strcmp(dirname, skip)))
        continue
    end
    
    if    ~strcmp( dirname,'.')           & ...
            ~strcmp( dirname,'..')        & ...
            ~strncmp( dirname,methodsep,1)& ...
            ~strcmp( dirname,'private')

        p = [p genpath2(fullfile(d,dirname), skips)]; % recursive calling of this function.
    end
end