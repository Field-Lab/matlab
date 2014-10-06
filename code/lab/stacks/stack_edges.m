function varargout = stack_edges(stack, islice, step)
% STACK_EDGES    Polygon coordinates representing bounds of stack image
% usage: [edges_xy] = stack_edges(stack, [islice, step])
% usage: [edges_x, edges_y] = stack_edges(stack, [islice, step])
%
% inputs:   islice      Defaults to 1
%           step        Step between point coordinates.  Defaults to 1.
%
% 2010-08 phli
%

if nargin < 2
    islice = 1;
end

if nargin < 3
    step = 1;
end


stack = load_slice_imfs(stack, islice);
imf = slice_imf(stack, islice);

edgelr = 1:step:imf.Width;
edgeud = 1:step:imf.Height;

edgesx = [edgelr                , imf.Width.*ones(1,length(edgeud)), fliplr(edgelr)                    , ones(1,length(edgeud)) 1];
edgesy = [ones(1,length(edgelr)), edgeud                           , imf.Height.*ones(1,length(edgelr)), fliplr(edgeud)         1];


if nargout == 1
    varargout{1} = [edgesx(:) edgesy(:)];
elseif nargout > 1
    varargout{1} = edgesx;
    varargout{2} = edgesy;
end