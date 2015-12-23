function [folders, names] = get_dataset_folder(ds_name, varargin)
% GET_DATASET_FOLDER follow the symbolic link on the server to find the folder
% containing the text files for a one dataset.
%
% the root one folder is assumed to be:
% 
%   '/snle/lab/Experiments/Array/Shared/one/'
%
% This file also contains the master list of all datasets that we are
% working with for the paper. it should be updated such that the list
% remains current. This file also verifies that all paths returned are
% valid. If any file/folder path is not valid, it returns an error.
% 
% [folder, names] = get_dataset_folder(name, varargin)
% 
% call get_dataset_folder('all')         to get all the folders (a cell array)
% call get_dataset_folder('plantain')    to get plantain's folder
%
% optional params, their default values, and what they specify:
%
% root_path         shown above      path to one text files
% suffix            none             add to end of paths (e.g. cones.txt)
%
% tamachado (6/24/09)

if ~isa(ds_name,'char')
   error('a dataset name (i.e. plantain) or all are the only valid arguments to this function');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

% specify list of optional parameters
p.addParamValue('root_path', '/snle/lab/Experiments/Array/Shared/one/');
p.addParamValue('suffix', '');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define all valid fruit names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%name = {'plantain','grapes','kiwi','mango','plum',...
%    'cherry','blueberry','butterfly','peach','apricot'};

name = {'plantain','grapes','kiwi',...
    'blueberry','peach','apricot'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save current working directory
current_directory = pwd;

% go to root path
cd(params.root_path);

% get indices of names we want
if strcmp(ds_name,'all')
    idr = 1:length(name);
else
    for ii = 1:length(name)
        if isequal(name{ii},ds_name), idr = ii; break; end
    end
end

% set up the return variables
names   = {name{idr}};
folders = cell(length(idr),1);
    
% follow link(s) to dataset(s) requested
for ii = 1:length(idr)
   % get the path
   cd([params.root_path '/' name{idr(ii)}]);
   folders{ii} = [pwd '/' params.suffix];
   
   % verify that it's a valid file/folder
   if ~exist(folders{ii}) %#ok<EXIST>
       cd(current_directory);
       error('%s does not exist!',folders{ii});
   end
   
end

% restore path to working directory
cd(current_directory);