function save_stacks(datarun, varargin)
% SAVE_STACKS       Save stack data to disk
% usage: save_stacks(datarun, opts)
%
% opts: strip   true    Whether to strip raw data before saving; in general
%                       if the path to the data is available, probably
%                       better to strip the data to avoid unnecessarily
%                       large save files
%
%       basepath    server_path()   Basepath to save; as of 2010-10 this defaults to marte Analysis
%
%       path        []      If set, this overrides path to save file to
%
%       file        stacks.mat  File name to save to
%
%       quiet       false
%
%       force       false   Overwrite existing data without prompting?
%
% Without options, this defaults to save to server_path/piece_name/stacks.mat
%
% 2010-09 phli
%


opts = inputParser;
opts.addParamValue('strip', true);
opts.addParamValue('basepath', server_path);
opts.addParamValue('path', []);
opts.addParamValue('file', 'stacks.mat');
opts.addParamValue('quiet', false);
opts.addParamValue('force', false);
opts.parse(varargin{:});
opts = opts.Results;


if ~isfield(datarun, 'stacks') || isempty(datarun.stacks)
    warning('No stack data to save');
    return
end

if ischar(opts.path)
    save_path = fullfile(opts.path);
else % Determine filename from basepath and piece name
    parsed_prefix = parse_rrs_prefix(datarun);
    if isempty(parsed_prefix)
        warning('Could not guess piece name');
        return
    end
    
    piece_name = parsed_prefix.piece_fullname;
    save_path = fullfile(opts.basepath, piece_name);
end
savepath = fullfile(save_path, opts.file);


if exist(savepath, 'file')
    overwrite = questdlg(['File ' savepath ' exists; overwrite?'], 'Overwrite existing stacks?', 'Okay', 'Cancel', 'Cancel');
    if strcmp(overwrite, 'Cancel')
        return
    end
end


stacks = datarun.stacks;
if opts.strip
    stacks = strip_stack_data(stacks);
end

if ~opts.quiet
    disp(['Saving stacks data to: ' savepath]);
end
save(savepath, 'stacks');