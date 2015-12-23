function datarun = init_stacks(datarun)
% INIT_STACKS   Initialize the stacks field if it does not yet exist
% usage: datarun = init_stacks(datarun)
%
% Initializes the stacks field with the array placeholder stack.  If 
% datarun.piece.array_id is set, this will also try to load in the array
% image stack.
%
% 2010-09 phli
%

% Do nothing if stacks is already initialized
if isfield(datarun, 'stacks') && ~isempty(datarun.stacks)
    return;
end

% See if we can get the array id from globals
if ~isfield(datarun, 'piece') || ~isfield(datarun.piece, 'array_id') || isempty(datarun.piece.array_id)
    datarun = load_globals(datarun);
    if isfield(datarun, 'globals') && ~isempty(datarun.globals) && datarun.globals.imageCalibrationParamsExists()
        icp = datarun.globals.getImageCalibrationParams();
        datarun.piece.array_id = icp.arrayID;
    end
end

% Initialize stacks; if we still don't have the array_id by now, punt
if isfield(datarun, 'piece') && isfield(datarun.piece, 'array_id') && ~isempty(datarun.piece.array_id)
    datarun.stacks = init_stacks(datarun.piece.array_id);
else
    datarun.stacks = init_stacks();
end