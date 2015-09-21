
function [x] = geninterp(f,T,timepts,samples)

1;
A = interpmtx(f,timepts(:),T);
x = A*samples(:);
%plot(timepts,samples,'o',1:T,x);
