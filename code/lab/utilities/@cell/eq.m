function bool = eq(cell1, cell2)
% CELL.EQ
%
% 2010-01-25 phli
%

warning('@cell.eq is deprecated; use ISEQUAL instead.');

if ~strcmp(class(cell1), class(cell2))
    bool = false;
    return;
end

if ~all(size(cell1) == size(cell2))
    bool = false;
    return;
end

for i = 1:numel(cell1)
    try
        sub1 = cell1{i};
        sub2 = cell2{i};
        
        if isempty(sub1) ~= isempty(sub2)
            bool = false;
            return;
        end
        
        if size(sub1) ~= size(sub2)
            bool = false;
            return;
        end
        
        if any(sub1(:) ~= sub2(:))
            bool = false;
            return;
        end
    catch e
        if strcmp(e.identifier, 'MATLAB:dimagree')
            bool = false;
            return;
        else
            throw(e);
        end
    end
end

bool = true;