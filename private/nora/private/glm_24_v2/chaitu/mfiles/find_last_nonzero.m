function cutoff = find_last_nonzero(x,tol)

x = x(:);
lastnonzero = find(abs(flipud(x)) > tol,1);
if (isempty(lastnonzero))
    cutoff = 1;
else
    cutoff = length(x) - lastnonzero;
end