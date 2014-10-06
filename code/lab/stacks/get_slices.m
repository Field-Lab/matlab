function slices = get_slices(stack, islices, reload)
% GET_SLICES    Get slices from image stack
% usage: slices = get_slices(stack, islices, reload)
%
% inputs:   stack     Image stack struct as described in Proposal.rtf
%           islices   Indices of slices to return
%           reload    Should the slices be freshly loaded even if they have already been cached?  Defaults to false.
%
% output: the image slices in cell array
%
% 2010-08 phli
%

if nargin < 3
    reload = false;
end

slices = cell(size(islices));
for i = 1:numel(islices)
    islice = islices(i);

    if ~reload && isfield(stack, 'data') && length(stack.data) >= islice && ~isempty(stack.data{islice})
        slices{i} = stack.data{islice};
        continue;
    end
    
    fullpath = slice_fullpath(stack, islice);
    page = stack.pages{islice};
    imf = slice_fullimf(stack, islice);
    if ~isempty(imf) && strcmp(imf(islice).Format, 'tif')     
        % Passing imf can make reading slices from stacked TIFFs much more efficient.
        % http://blogs.mathworks.com/steve/2009/04/02/matlab-r2009a-imread-and-multipage-tiffs/
        % ToDo: Except I'm not really seeing these efficiency wins for stacks of <= 4000...?
        im = imread(fullpath, page, 'Info', imf);
    elseif page > 1
        im = imread(fullpath, page);
    else
        im = imread(fullpath);
    end
    
    slices{i} = im;
end