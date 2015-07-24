% Function that permutes the entries of a cell array

function C2 = permute_cellarray(C,neworder)

C2 = cell(size(C));

for j=1:length(C)
    C2{j} = C{neworder(j)};
end
