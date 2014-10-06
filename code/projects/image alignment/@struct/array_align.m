function datarun = array_align(datarun, alive_index, varargin)
% ARRAY_ALIGN
% usage: datarun = array_align(datarun, alive_index, opts)
%
% Align an image to array coordinates, using simpler corners method.  For
% this to work, the image has to include the array corners and not be
% warped.
%
% See also align_alive_to_array cell/align_alive_to_array
%
% 2010-08 phli
%

datarun.stacks = align_alive_to_array(datarun.stacks, alive_index, datarun, varargin{:});