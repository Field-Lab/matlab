function datarun = stack_align(datarun, varargin)
% STACK_ALIGN
% usage: datarun = stack_align(datarun, input_index, base_index, opts)
%
% See also cell/stack_align
%
% 2010-08 phli
%

datarun.stacks = stack_align(datarun.stacks, varargin{:});