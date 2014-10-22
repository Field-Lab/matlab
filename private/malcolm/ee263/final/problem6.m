% problem 6
clear all; close all; clc;

load('facial_recognition_data.mat');

figure();
colormap gray;
imagesc(Y);
axis off;

% 6b: find Fk that is closest to Y
E_f = nan(N,1);
for k = 1:N
    Fk = F{k};
    E_f(k) = norm(Y-Fk,'fro')^2;
end
k = find(E_f == min(E_f));
figure();
colormap gray;
imagesc(F{k});
axis off;
% k =
%     59

% 6c
f = nan(m*n,N);
for k = 1:N
    f(:,k) = F{k}(:);
end
fbar = (1/N)*sum(f,2);
X = f - repmat(fbar,1,N);
[U,S2] = eigs(X*X',r);
p = nan(r,1);
sigma2 = diag(S2);
sigma_sum = trace(X*X');
for i = 1:r
    p(i) = sum(sigma2(1:i))/sigma_sum;
end
figure();
stem(p);
xlabel('i');
ylabel('p_i');
title('Problem 6c');
% p(r)
%     0.8348

% 6d
nrow = 4;
ncol = 5;
figure();
colormap gray;
for k = 1:r
    subplot(nrow,ncol,k);
    imagesc(reshape(U(:,k),m,n));
    axis off;
end

% 6e
y = Y(:);
figure();
colormap gray;
for d = 1:r
    subplot(nrow,ncol,d);
    yhatd = reshape(fbar+U(:,1:d)*U(:,1:d)'*(y-fbar),m,n);
    imagesc(yhatd);
    axis off;
end

% 6f
Ur = U(:,1:r);
E_u = nan(N,1);
for k = 1:N
    E_u(k) = norm(Ur'*(f(:,k)-y));
end
k = find(E_u == min(E_u));
figure();
colormap gray;
imagesc(F{k});
axis off;
% k = 
%     91