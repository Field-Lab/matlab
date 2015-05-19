% problem 5
close all; clc;

sigma = svd(votes);
figure();
stem(sigma);
xlabel('i');
ylabel('\sigma_i(votes)');
title('Problem 5a');
f2 = sum(sigma(1:2).^2)/sum(sigma.^2);
% f2 =
%     0.7307

tol = 1e-3;
A = votes;
Ahat = A;
first = true;
while (norm(A-Ahat) > tol) || first
    first = false;
    A = Ahat;
    [U,S,V] = svd(A);
    Ahat = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
    Ahat(find(votes==0)) = sign(Ahat(find(votes==0)));
    Ahat(find(votes~=0)) = votes(votes~=0);
end
% sum(sum(A==1))
%     44934
% sum(sum(A==-1))
%     19632

% U, S, and V from above are the svd of A
u1 = U(:,1);
u2 = U(:,2);
z1 = nan(102,1);
for i = 1:102
    z1(i) = 2*(senators{i,3} == 'D') + (senators{i,3} == 'I');
end
spatial_plot(u1, u2, z1, 3, eye(3));
xlabel('u1');
ylabel('u2');
title('Political party');

z2 = nan(102,1);
for i = 1:102
    z2(i) = sum(A(i,:) == sign(sum(A)));
end
spatial_plot(u1, u2, z2, 10);
xlabel('u1');
ylabel('u2');
title('Rate of agreement with majority');

v1 = V(:,1);
v2 = V(:,2);
z3 = sum(A);
spatial_plot(v1, v2, z3, 10);
xlabel('v1');
ylabel('v2');
title('Total support');

z4 = nan(633,1);
R = sum(z1 == 0);
D = sum(z1 == 2);
for j = 1:633
    z4(j) = (1/R)*sum(A(z1==0,j)) - (1/D)*sum(A(z1==2,j));
end
spatial_plot(v1, v2, z4, 10);
xlabel('v1');
ylabel('v2');
title('Partisan support');