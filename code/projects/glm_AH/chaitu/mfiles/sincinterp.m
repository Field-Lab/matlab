function [x timepts] = sincinterp(taxis,timepts,samples)

A = sincinterpmtx(taxis(:),timepts(:));
x = A*samples;
if 0
    plot(timepts,samples,'o',taxis,x,'k.-');
end
