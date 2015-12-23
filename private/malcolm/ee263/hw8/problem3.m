% problem 3
close all; clc

% construct the node-incidence matrix
B = zeros(m,n);
counter = 1;
for i = 1:n
    for j = i:n
        if A(i,j)
            B(counter,min(i,j)) = 1;
            B(counter,max(i,j)) = -1;
            counter = counter+1;
        end
    end
end

[~,~,v1] = svds(B,1);
x = sign(sqrt(n)*v1);

% calculate number of cuts using the SVD method
cuts_svd = norm(B*x)^2

% generate random partitions
N = 1000;
cuts_rand = nan(N,1);
for i = 1:N
    xrand = sign(rand(n,1)-0.5);
    cuts_rand(i) = norm(B*xrand)^2;
end
max_cuts_rand = max(cuts_rand)

% cuts_svd =
%    268
% 
% max_cuts_rand =
%    196