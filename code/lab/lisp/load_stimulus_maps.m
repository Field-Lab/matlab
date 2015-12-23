function stimstruct = load_stimulus_maps(stimstruct, piece, varargin)
% usage: stimstruct = load_stimulus_maps(stimstruct)
% STIMSTRUCT.SPEC.INDEX_MAP is generated from stim lisp output, so multiple 
% formats must be handled.
%
% See also READ_STIM_LISP_OUTPUT
%
% 2012-06, phli
%

opts = inputParser();
opts.addParamValue('parsepoly', true);
opts.parse(varargin{:});
opts = opts.Results;

if ~isfield(stimstruct.spec, 'index_map'), return; end

mappath = stimstruct.spec.index_map;
mappath = strrep(mappath, ':', filesep());


% Mappath often comes out of stim lisp output with a filesep() on the front
% even though it's not a root directory, so can't use that as check for
% root directory.


% First check in the stimulus_maps_path
stimmappath = add_base_path(mappath, stimulus_maps_path());
if exist(stimmappath, 'file')
    stimstruct = load_stimulus_maps_(stimstruct, stimmappath, opts);
    return
end


% Now check analysis directories
% Sometimes stim lisp output puts the piece on the front of the mappath,
% but this doesn't work nicely with the mappath format used in the analysis
% directories, so try to pull this off the front.
rgxp = ['\' filesep '(\d{4}-\d{2}-\d{2}-\d)' '(\' filesep '.*)' ];
tokens = regexp(mappath, rgxp, 'tokens');
if length(tokens) == 1 % Have to check what it means if > 1...
    tokens = tokens{1};
    
    if isempty(piece)
        piece = tokens{1};
    end
    
    if strcmp(piece, tokens{1})
        mappath = tokens{2};
    end
end
analmappath = fullfile(server_path, piece, 'stimuli', mappath);
if exist(analmappath, 'file')
    stimstruct = load_stimulus_maps_(stimstruct, analmappath, opts);
    return
else
    fprintf('Can''t find %s\n', analmappath);
end



function stimstruct = load_stimulus_maps_(stimstruct, mappath, opts)
stimstruct.mappath = mappath;

if exist(mappath, 'dir')
    % It is a directory of map files
    mapfiles = dir([mappath '/map-*.txt']);
    for i = 1:length(mapfiles)
        stimstruct.mapims{i} = sparse(dlmread([mappath '/' mapfiles(i).name]));
    end
else
    % It is a single filename
    stimstruct.mapims{1} = sparse(dlmread(mappath));
end

for i = 1:length(stimstruct.mapims)
    stimstruct.mapnyc{i} = map2manhattan(stimstruct.mapims{i}, 'quiet', true);
    if ~opts.parsepoly, continue; end;
    
    for j = 1:length(stimstruct.mapnyc{i})
        stimstruct.mapnycpoly{i}{j} = segs2poly(stimstruct.mapnyc{i}{j});
    end
end