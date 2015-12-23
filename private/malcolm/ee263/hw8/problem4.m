% problem 4b
close all; clc;

[dinv,~,a] = svds(G,1);
% renormalize
% check that all entries of dinv are nonzero
if all(dinv)
    Ghat = dinv*a';
    d = 1./dinv;
    c = mean(d);
    d = d/c
    a = a/c;
    Jopt = sqrt(1/(m*n))*norm(G-Ghat,'fro')
    ratio = Jopt/norm(G,'fro')
else
    disp('some entries of dinv are 0, so cannot invert');
end

% d =
%     0.9429
%     1.2780
%     0.9015
%     0.9197
%     0.7729
%     1.0418
%     1.1433
%
% Jopt =
%    15.8400
% 
% ratio =
%     0.0427
