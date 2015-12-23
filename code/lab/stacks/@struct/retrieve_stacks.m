function loaded = retrieve_stacks(datarun, varargin)
% RETRIEVE_STACKS       Load stack data for given datarun from disk
% usage: loaded = retrieve_stacks(datarun, opts)
%
% opts: basepath    server_path()       As of 2010-10 this defaults to marte Analysis
%       path        []                  If set, overrides path to stack saved data file
%       file        'stacks.mat'        Filename of saved stack data
%
% Without any options, defaults to load from server_path/piece_name/stacks.mat
%
% 2010-10 phli
%

opts = inputParser;
opts.addParamValue('basepath', server_path);
opts.addParamValue('path', []);
opts.addParamValue('file', 'stacks.mat');
opts.parse(varargin{:})
opts = opts.Results;

if ischar(opts.path)
    load_path = opts.path;
else % Guess loadpath from basepath and piece name
    parsed_prefix = parse_rrs_prefix(datarun);
    if isempty(parsed_prefix)
        error('Could not guess piece name');
    end
    
    piece_name = parsed_prefix.piece_fullname;
    load_path = fullfile(opts.basepath, piece_name);
    if ~isdir(load_path)
        warning(['Not a valid path: ' load_path]);
        return
    end
end

% Look for stacks file, add file name to dir load_path
files = dir(load_path);
files = {files.name};
if ~ismember(opts.file, files)
    warning(['No saved stacks file found: ' load_path '/' opts.file]);
    return
end
load_path = fullfile(load_path, opts.file);

% Load stacks
loaded = load(load_path);