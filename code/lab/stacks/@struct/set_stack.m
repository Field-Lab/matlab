function datarun = set_stack(datarun, stack, varargin)
% SET_STACK   Assign stack to datarun
% usage: datarun = set_stack(datarun, stack, stack_spec)
%        datarun = set_stack(datarun, stack)
%
% Stack spec is as described in PARSE_STACK_INDEX.  If the stack has a
% "name" field that matches a named stack index in STACK_LABELS, then that
% value can be used and the stack_spec argument can be left off.
%
% See also: PARSE_STACK_INDEX, STACK_LABELS
%
% 2010-09 phli
%

datarun = init_stacks(datarun);
datarun.stacks = set_stack(datarun.stacks, stack, varargin{:});