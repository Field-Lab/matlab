function bool = stackempty(datarun, varargin)
% STACK_EMPTY   Return boolean indicating if given stack is empty
% usage: bool = stackempty(datarun, stack_spec)
% e.g.:  bool = stackempty(datarun, {2,2});
%        bool = stackempty(datarun, 2, 2);
%        bool = stackempty(datarun, 'alive_montage_hires');
%
% Returns true (i.e. empty) if stacks does not have high enough dimension
% to accomodate stack_spec.
%
% 2010-10 phli
%
bool = true;

if ~isfield(datarun, 'stacks')
    return;
end

bool = stackempty(datarun.stacks, varargin{:});