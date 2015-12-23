function M = cellnest2matrix(N)
% Takes nested cell vectors and makes them into multidimensional cell
% matrices.  For example:
%
% cellnest2matrix({{1 2 3} {4 5 6} {7 8 9}})
% 
% ans = 
% 
%     [1]    [2]    [3]
%     [4]    [5]    [6]
%     [7]    [8]    [9]
%
% But this will work for any number of levels of nesting.  Of course, if
% the nesting structure doesn't conform to nests of cell vectors, then
% things will go wrong.  Can handle slightly more flexible situation of one
% vector in a particular level being different length than another.
%
% Being able to do this simplifies writing recursive functions that build
% up multidimensional matrices; instead of worrying about the
% dimensionalities when writing the recursion, you can just return nested
% cell arrays, which is simpler, and then apply this at the end.
    
if islastlevel(N),
    M = N;
    return
end

% Kludge to handle 2D representation of vectors at lowest level
ispenultimatelevel = any(cellfun(@islastlevel, N));

maxsize = [];
maxlen = length(maxsize);
for i = 1:length(N)
    temp{i} = cellnest2matrix(N{i});
    tempsize = size(temp{i});
    if ispenultimatelevel
        tempsize = squeeze2D(tempsize);
    end
    
    for j = 1:length(tempsize)
        if j > maxlen || tempsize(j) > maxsize(j)
            maxsize(j) = tempsize(j);
        end
    end
    maxlen = length(maxsize);
end

M = cell([length(N) maxsize]);
S.type = '()';
for i = 1:length(temp)
    tempsize = size(temp{i});
    if ispenultimatelevel
        tempsize = squeeze2D(tempsize);
    end
    
    S.subs = {i};
    for j = 1:length(tempsize)
        S.subs{j+1} = 1:tempsize(j);
    end
    
    M = subsasgn(M, S, temp{i});
end


function bool = islastlevel(N)
bool = false;
for i = 1:length(N)
    if ~iscell(N{i})
        bool = true;
        return
    end
    
    if length(N{i}) ~= numel(N{i})
        bool = true;
        return
    end
end


% Handle special case of size vector for 2D; don't add extra dimension
function s = squeeze2D(s)
if length(s == 2) && s(1) == 1
    s = s(2);
end