function [x y c] = stack_xy(stack)
% STACK_XY
% usage: [x y c] = stack_xy(stack)
%
% Only tested to work on TIFF based stacks; that's mostly what we use; add
% more formats as they come up.
%
% 2010-12 phli
%

x = 0;
y = 0;
c = 0;

if isempty(stack)
    return
end


imf = slice_imf(stack);
if ~isempty(imf)
    x = imf.Width;
    y = imf.Height;
    c = imf.SamplesPerPixel;
    return
end


slice = get_slice(stack, 1);
[y x c] = size(slice);