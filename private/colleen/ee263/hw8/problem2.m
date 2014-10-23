% problem 2
close all; clc

% 2a
Atilde = nan(m,n);
for i = 1:n
    Atilde(:,i) = A(:,i)/norm(A(:,i));
end
[u,s,v] = svd(Atilde);
figure();
stem(diag(s));
title('Problem 2a');
xlabel('Singular values of Atilde');

% 2b
% query for "students" (i=53)
qtilde = zeros(m,1);
qtilde(53) = 1;
c = Atilde'*qtilde;
[~,idx] = sort(c,'descend');

% 2c
% low rank search results
idx_approx = nan(n,4);
for i = 1:4
    r = 32/(2^(i-1));
    Aapprox = u(:,1:r)*s(1:r,1:r)*v(:,1:r)';
    c_approx = Aapprox'*qtilde;
    [~,idx_approx(:,i)] = sort(c_approx,'descend');
end

% show top 5 results for all searches
% >> [idx(1:5) idx_approx(1:5,:)]
%
% ans =
%
%    106   106   106   115   115
%    105   105   107   106   105
%    107   107   105   120   107
%    115   115   115   107    66
%    111   111   111   111    63