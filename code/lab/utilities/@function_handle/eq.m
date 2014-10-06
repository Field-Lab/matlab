function are_equal = eq(fh1, fh2)
% FUNCTION_HANDLE.EQ
%
% 2010-05 phli
%

if ~strcmp(class(fh1), class(fh2))
    are_equal = false;
    return;
end

fs1 = func2str(fh1);
fs2 = func2str(fh2);
are_equal = strcmp(fs1, fs2);