function [e1,e2]=ext(d)
%extreme

[t, e2]=max(abs(d));
e1=d(e2);