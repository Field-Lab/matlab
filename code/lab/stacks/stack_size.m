function siz = stack_size(stack)
% STACK_SIZE
% usage: [x y c l] = stack_size(stack)
%
%
% 2010-12 phli
%

[x y c] = stack_xy(stack);
l = stack_length(stack);
siz = [y x c l];