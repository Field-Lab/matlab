function bool = eq(struct1, struct2)
% STRUCT.EQ
%
% 2010-01-25 phli

if class(struct1) ~= class(struct2)
    bool = false;
    return;
end

% Structs can be arrays
if ~all(size(struct1) == size(struct2))
    bool = false;
    return;
end

for i = 1:numel(struct1)
    fieldnames1 = fieldnames(struct1(i));
    if sort(fieldnames1) ~= sort(fieldnames(struct2(i)))
        bool = false;
        return;
    end
    
    try
        for j = 1:length(fieldnames1)
            sub1 = struct1(i).(fieldnames1{j});
            sub2 = struct2(i).(fieldnames1{j});
            if isempty(sub1) ~= isempty(sub2)
                bool = false;
                return;
            end
            
            if any(sub1 ~= sub2)
                bool = false;
                return;
            end
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