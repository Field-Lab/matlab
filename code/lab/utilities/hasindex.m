function bool = hasindex(A, indices)
% HASINDEX  Check that array A is large enough to have the given index/indices
% usage: hasindex(A, indices)

% Just convert array of indices to cell array
indices = num2cell(indices);

bool = true;
try
    sub2ind(size(A), indices{:});
catch e
    if strcmp(e.identifier, 'MATLAB:sub2ind:IndexOutOfRange');
        bool = false;
    else
        rethrow(e);
    end
end