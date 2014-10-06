function typestr = stack_datatype(stack)
% STACK_DATATYPE
% usage: typestr = stack_datatype(stack)
%
% Not very efficient
%
% 2010-12 phli
%

typestr = class(get_slice(stack));