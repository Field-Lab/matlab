% problem 4
close all; clc;

% 4a
y1true = A1*xtrue;
y2true = A2*xtrue;
xhat = A1\y1tilde;
y2hat = A2*xhat;
ex = norm(xtrue-xhat)/norm(xtrue);
ey2 = norm(y2true-y2hat)/norm(y2true);
% ex =
%    1.6944e+04
% ey2 =
%     0.2443

% 4b
A = [A1; A2];
[U,S,V] = svd(A);
figure();
stem(diag(S));
xlabel('i');
ylabel('\sigma_i(A)');
title('Problem 4b');

% 4c
p = 3;
U11 = U(1:m1,1:p);
U21 = U(m1+1:m1+m2,1:p);
S1 = S(1:p,1:p);
V1 = V(1:n,1:p);
ztrue = V1'*xtrue;
zhat = (U11*S1)\y1tilde;
ez = norm(ztrue-zhat)/norm(ztrue);

% zhat =
%     0.0266
%     0.7281
%     0.8991
% ztrue =
%     0.0284
%     0.7333
%     0.8998
% ez =
%     0.0047

% 4d
y2hatz = U21*S1*zhat;
ey2z = norm(y2true-y2hatz)/norm(y2true);
% ey2z =
%     0.0052