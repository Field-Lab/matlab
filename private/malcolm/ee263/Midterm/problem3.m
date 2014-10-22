% problem 3a
m = length(x);
X = [x ones(m, 1)];
a_ls = X\y;
plot(x,y,'.','MarkerSize',20)
h1=refline(a_ls(1),a_ls(2));
set(h1,'LineStyle','--')
xlabel('x')
ylabel('y')
title('Problem 3a')

% a_ls =
%     1.8385
%    15.9447

% problem 3b
% iteratively reweighted least squares
a_mpe_new = a_ls;
tol = 1e-3;
a_mpe = [0;0];
counter = 0;
while norm(a_mpe_new-a_mpe) > tol
    a_mpe = a_mpe_new;
    W = sqrt(diag(weights(x,y,a_mpe)));
    a_mpe_new = (W*X)\(W*y);
    counter = counter+1;
end
h2=refline(a_mpe(1), a_mpe(2));
legend({'data','LS fit','MPE fit'},'location','northwest')

% a_mpe =
%     1.9700
%    10.6769

% problem 1c
a = [2;8];
rmse_ls = sqrt(sum((a-a_ls).^2)/2);
rmse_mpe = sqrt(sum((a-a_mpe).^2)/2);

% rmse_ls =
%     5.6189  
% rmse_mpe =
%     1.8930
    