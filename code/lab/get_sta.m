function [sta,sta_stored] = get_sta(datarun,cell_id,varargin)
% get_sta     return the STA for a cell.  If needed, read the STA from disk.
%
% usage:  [sta,sta_stored] = get_sta(datarun,cell_id,varargin)
%
% arguments:  datarun - datarun struct
%             cell_id - cel id
%            varargin - struct or list of optional parameters (see below)
%
% outputs:        sta - 4D matrix (y,x,color,time)
%          sta_stored - boolean.  if true, the sta was located in datarun.  if false, it needed to be loaded.
%
%
% optional params, their default values, and what they specify:
%
% frames        ':'     which frames to load
%
%
% gauthier 2008-10
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('frames', ':');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION

% get cell index
cell_index = get_cell_indices(datarun,cell_id);

% if the STA is loaded, use it
if isfield(datarun,'stas') && isfield(datarun.stas,'stas') && length(datarun.stas.stas) >= cell_index && ~isempty(datarun.stas.stas{cell_index})
    sta = datarun.stas.stas{cell_index};
    sta_stored = true;
else
    % otherwise, load the STA

    if ~isfield(datarun, 'stimulus') || ~isfield(datarun.stimulus, 'independent') || ~ischar(datarun.stimulus.independent)
        datarun = load_sta(datarun, 'load_sta', []);
    end
    
    % if the java sta object exists, use it
    if isfield(datarun,'stas') && isfield(datarun.stas,'java_sta') && ~isempty(datarun.stas.java_sta)
        java_sta = datarun.stas.java_sta;
    else
        % otherwise, load the java sta object from disk
        java_sta = load_java_sta(datarun,'verbose',false);
    end

    % get the STA
    sta = sta_from_java_stas(java_sta, cell_id, params.frames, datarun.stimulus.independent);
    
    sta_stored = false;
end
