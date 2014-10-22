% problem 1
close all; clc;

% define the matrix A
A = nan(n,n);
for i = 1:n
    for j = 1:n
        if i==j
            A(i,j) = 1-sum(K(i,1:n ~=i));
        else
            A(i,j) = K(j,i);
        end
    end
end

% diagonalize A
[V,D] = eig(A);
lambda = diag(D);
if det(V) == 0
    fprintf('A is not diagonalizable\n');
    return;
else
    fprintf('A is diagonalizable\n');
end
W = inv(V);

if sum(abs(lambda)==max(abs(lambda))) ~= 1
    disp('There is not a uniquely dominant eigenvalue.');
    return;
else
    disp('There is a uniquely dominant eigenvalue.');
end
v1 = V(:,abs(lambda)==max(abs(lambda)));
w1 = W(abs(lambda)==max(abs(lambda)),:)';

% compute r
r = sort(abs(lambda),'descend');
r = r(2);

% simulate for t = 0:500
T = 500;
x0 = eye(n,n);
x0 = x0(:,7);
xbar = v1*w1'*x0;
x = nan(n,T+1);
x(:,1) = x0;
e = nan(T+1,1);
e(1) = norm(x0-xbar);
for t = 1:T
    x(:,t+1) = A*x(:,t);
    e(t+1) = norm(x(t+1)-xbar);
end

% sanity check: sum(xbar) = 1
% sum(xbar)
%     1

% make plots
figure();
subplot(2,1,1);
plot(0:T,x(10,:));
hline = refline(0,xbar(10));
set(hline,'LineStyle','--');
ylim([0 0.12]);
ylabel('x_{10}(t)');
title('Problem 1d');
subplot(2,1,2);
semilogy(0:T,e);
xlabel('t');
ylabel('||e(t)||');