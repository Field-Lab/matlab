function stack = load_slices(stack, islices, reload)
% LOAD_SLICES
% usage: stack = load_slices(stack, islices, reload)
%
% inputs:   stack       Image stack struct as described in Proposal.rtf
%           islices     Indices of slices to load
%           reload      Whether to reload if already cached; default false
%
% For retrieving, see also get_slices get_slice
%
% 2010-08 phli
%

if nargin < 2
    islices = 1:stack_length(stack);
end

if nargin < 3
    reload = false;
end

for islice = islices
    if isfield(stack, 'data') && length(stack.data) >= islice && ~isempty(stack.data{islice}) && ~reload
        continue;
    end
    
    im = get_slices(stack, islice, reload);
    stack.data(islice) = im;
end