function flat = flatten(in)
% FLATTEN   Flatten cell array
%
% 2012-06, phli
%

if ~iscell(in)
    flat = in;
    return
end


flat = {};
for i = 1:numel(in)
    flat = [flat flatten(in{i})];
end