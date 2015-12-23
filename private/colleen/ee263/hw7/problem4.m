% problem 4b
clc; close all;

% b is just the mean of Y
b = mean(Y,2);

% now use the svd to calculate Ax
Yshift = Y-repmat(b,1,100);
[U,S,V] = svd(Yshift);

% try many different values of m
m = 1:10;
J = nan(length(m),1);
for i = m
    A = U(:,1:i)*S(1:i,1:i);
    X = V(:,1:i)';
    J(i) = sqrt((1/N)*norm(Yshift-A*X,'fro')^2);
end

% plot J vs m
figure();
plot(m,J);
xlabel('m');
ylabel('J');
title('Problem 4b-1');

% Plot shows that m = 3.
%
% A =
%    -2.2978   -0.3194   -0.0166
%    -2.9712   -0.4408   -0.4928
%    -1.3005    0.4375   -0.1360
%    -2.7130    0.6446   -0.0780
%    -0.7816    0.8584    0.0727
%    -1.5121   -0.2081    0.4146
%    -1.6140    0.1006    0.5250
%    -1.8012    0.8874   -0.0420
%    -1.3734   -0.6178   -0.3576
%    -2.0961   -0.7360    0.4420
%
% b =
%     0.0629
%     0.0811
%     0.0511
%     0.1005
%     0.0457
%     0.0411
%     0.0525
%     0.0788
%     0.0273
%     0.0460

% make stem plot
A = U(:,1:3)*S(1:3,1:3);
X = V(:,1:3)';
norm_error = sqrt(sum((Yshift-A*X).^2));
figure();
stem(sort(norm_error,'descend'));
ylabel('norm of error');
xlabel('index');
title('Problem 4b-2');

