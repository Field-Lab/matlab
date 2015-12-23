% problem 5b
close all; clc;

% initialize the estimate Ahat
Ahat = ones(m,n)*mean(Aknown);
for i = 1:length(Aknown)
    Ahat(K(i,1),K(i,2)) = Aknown(i);
end

% perform iteration kmax times
kmax = 300;
norm_error = nan(kmax,1);
for i = 1:kmax
    [U,S,V] = svd(Ahat);
    Ahat = U(:,1:r)*S(1:r,1:r)*(V(:,1:r)');
    for j = 1:length(Aknown)
        Ahat(K(j,1),K(j,2)) = Aknown(j);
    end
    norm_error(i) = norm(A-Ahat,'fro');
end

% plot norm error vs iteration
figure();
plot(1:kmax,norm_error);
xlabel('iteration');
ylabel('norm error');
title('Problem 5b');