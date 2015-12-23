% problem 3
close all; clc;

Xbar = mean(X,2);
Xtilde = X - repmat(Xbar,1,N);
[U,~,~] = svd(Xtilde);
Q = U(:,1:2);
Xhat = Q'*Xtilde;

VRn = (1/N)*norm(Xtilde,'fro')^2;
VS = (1/N)*norm(Q*Xhat,'fro')^2;
VStilde = (1/N)*norm(Xtilde(1:2,:),'fro')^2;

% VS/VRn
%     0.3006
% VStilde/VRn
%     0.0142

subplot(2,1,1);
Xproj = Q'*X;
plot(Xproj(1,:),Xproj(2,:),'x');
xlabel('q1');
ylabel('q2');
title('Problem 3b');
subplot(2,1,2);
plot(X(1,:),X(2,:),'x');
xlabel('e1');
ylabel('e2');