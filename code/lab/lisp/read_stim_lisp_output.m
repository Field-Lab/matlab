function parsed = read_stim_lisp_output(path_or_piece, blank_or_stim, mappath, parsepoly)
% READ_STIM_LISP_OUTPUT     Attempt to read in and parse stimulus Lisp output
% usage: parsed = read_stim_lisp_output(filepath, path_or_piece, blank_or_stim, mappath, parsepoly)
%
% These usually live in files like "s01".  There are different stimulus
% types; only a few of them have parsers set up.  This will attempt to
% parse as far as it knows how.
%
% If filepath is not absolute, it will be assumed to live under
% SERVER_DATA_PATH.
%
% See also: SERVER_DATA_PATH, LISP2CA
%
% 2011-05 phli
% This is redundant with some of Martin's functions.
%

if nargin < 2 || isempty(blank_or_stim)
    piece = '';
    stim = '';
    filepath = path_or_piece;
else
    piece = path_or_piece;
    stim = blank_or_stim;
    filepath = fullfile(piece, 'stimuli', stim);
end
if nargin < 3, mappath = []; end
if nargin < 4, parsepoly = true; end


% Add SERVER_DATA_PATH unless already absolute path
filepath = add_base_path(filepath, server_path);

% File exists?
if ~exist(filepath, 'file')
    error('Could not find stimulus file %s.', filepath);
end

% Read lines from file
lines = {};
f = fopen(filepath, 'r');
tline = fgetl(f);
while ischar(tline) && ~isempty(tline)
    lines{end+1} = tline;
    tline = fgetl(f);
end
fclose(f);

% Convert Lisp lines into cell arrays, adding them to a master cell array
cells = {};
for i = 1:length(lines)
    line = lines{i};
    
    % There is one complication: often multiple () () () of Lisp are on a
    % single line.  When converted to {} {} {} this will not evaluate as
    % valid Matlab.  So if LISP2CA fails, we assume that the problem is
    % that the line needs to be wrapped in an extra outer ().  We then cut
    % this added {} back out of the entry in the master cell array.
    try
        cells{end+1} = lisp2ca(line);
    catch
        tempcells = lisp2ca(['(' line ')']);
        cells(length(cells)+1:length(cells)+length(tempcells)) = tempcells;
    end
end


% The first sexp of the stimulus output Lisp should be a general
% stimulus information sexp
spec = cells{1};
if ~strcmpi(spec{1}, 'type')
    warning('Unrecognized stimulus Lisp output format; first symbol token should be "TYPE"; returning raw cell array');
    parsed = cells;
    return
end


% The stimulus specification sexp should be in a paired key value format:
% (:type TYPE key1 value1 key2 value2 etc.)
parsed.spec = parse_stim_lisp_output_kv(spec);
rest = cells(2:end);

% Parse the rest of the sexps specifically according to type
switch lower(parsed.spec.type)
    case 'pulse'
        parsed = parse_pulse_output(parsed, rest);
    case 'pulse-sequence'
        parsed = parse_pulse_sequence_output(parsed, rest);
    case 'movie'
        parsed = parse_movie_output(parsed, rest);
    case 'drifting-sinusoid'
        parsed = parse_drifting_sinusoid_output(parsed, rest);
    case 'reversing-sinusoid'
        parsed = parse_reversing_sinusoid_output(parsed, rest);
    otherwise
        warning(['Unrecognized stimulus type ' parsed.spec.type '; returning raw cell arrays']);
        parsed.cells = cells(2:end);
end


% Attempt to load stimulus map if any
if ~isempty(mappath), parsed.spec.index_map = mappath; end
parsed = load_stimulus_maps(parsed, piece, 'parsepoly', parsepoly);


function strct = parse_stim_lisp_output_kv(ca)
for i = 1:2:length(ca)
    fldnam = lower(ca{i});
    fldnam = strrep(fldnam, '-', '_');
    strct.(fldnam) = ca{i+1};
end


% Generic sinusoid
function parsed = parse_sinusoid_output(parsed, cells)
for i = 1:length(cells)
    parsed.sinusoids(i) = parse_stim_lisp_output_kv(cells{i});
end
parsed.spatialperiods  = [parsed.sinusoids.spatial_period];
parsed.temporalperiods = [parsed.sinusoids.temporal_period];

% Below is simpler, but for parse_stim_rgbs / numcellunique we need this
parsed.rgbs = arrayfun(@(S)([S.rgb{:}]), parsed.sinusoids, 'UniformOutput', false);
%parsed.rgbs = cell2mat(vertcat(parsed.sinusoids.rgb));


% DRIFTING_SINUSOID
function parsed = parse_drifting_sinusoid_output(parsed, cells)
parsed = parse_sinusoid_output(parsed, cells);
parsed.directions = [parsed.sinusoids.direction];


% REVERSING_SINUSOID
function parsed = parse_reversing_sinusoid_output(parsed, cells)
parsed = parse_sinusoid_output(parsed, cells);
parsed.spatialphases = [parsed.sinusoids.spatial_phase];


%PULSE
function parsed = parse_pulse_output(parsed, cells)
for i = 1:length(cells)
    parsed.pulses(i) = parse_stim_lisp_output_kv(cells{i});
end

firstpulse = parsed.pulses(1);
if isfield(firstpulse, 'rgbs')
    parsed.numcones = length(firstpulse.rgbs);
    parsed.rgbs = collect(parsed.pulses, @(s) (cell2mat(collect(s.rgbs', @cell2mat))));
elseif isfield(firstpulse, 'rgb')
    parsed.numcones = 1;
    parsed.rgbs = arrayfun(@(strct)(cell2mat(strct.rgb)), parsed.pulses, 'UniformOutput', false);
end

if isfield(parsed.pulses, 'map')
    parsed.maps = [parsed.pulses.map];
    parsed.mapindices = unique(parsed.maps);
end


% PULSE_SEQUENCE
function parsed = parse_pulse_sequence_output(parsed, cells)
for i = 1:length(cells)
    parsed.pulses(i) = parse_stim_lisp_output_kv(cells{i});
    parsed.pulses(i).rgb = cell2mat(parsed.pulses(i).rgb);
end
parsed.maps = [parsed.pulses.map];
parsed.rgbs = vertcat(parsed.pulses.rgb);
parsed.mapindices = unique(parsed.maps);


% MOVIE
function parsed = parse_movie_output(parsed, cells)
for i = 1:length(cells)
    kv(i) = parse_stim_lisp_output_kv(cells{i});
end
parsed.start_frames = [kv.start_frame];
