function out = structwhos(S)
% STRUCTWHOS    Run whos on the internal fields of struct
% usage: out = structwhos(S)
%
% Helpful for comparing two versions of a struct for size, dataformat, etc.
% Like WHOS, if there's no argout it will print to screen, otherwise will
% return results to OUT silently.
%
% 2012-07-18, phli
%

fields = fieldnames(S);
for f = 1:length(fields)
    eval(sprintf('%s = S.(''%s'');', fields{f}, fields{f}));
end

if nargout == 0
    whos(fields{:});
else
    out = whos(fields{:});
end