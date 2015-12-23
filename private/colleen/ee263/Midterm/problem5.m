% problem 5b
% define Aold_tilde, Anew_tilde, and etilde
etilde = reshape(eye(n),n*n,1);
Aold_tilde = nan(n*n,m*n);
Anew_tilde = nan(n*n,m*n);
for i = 1:n
    Aold_tilde((i-1)*n+1:i*n,:) = [zeros(n,(i-1)*m) Aold' zeros(n,m*(n-i))];
    Anew_tilde((i-1)*n+1:i*n,:) = [zeros(n,(i-1)*m) Anew' zeros(n,m*(n-i))];
end
C = [Aold_tilde'*Aold_tilde Anew_tilde'; 
    Anew_tilde zeros(n*n)];
Btilde_lambda = C\[Aold_tilde'*etilde; etilde];
Btilde = Btilde_lambda(1:m*n);
B = reshape(Btilde,m,n)';
% calculate J using this B and the pseudoinverse of Anew
Jopt = norm(B*Aold-eye(n),'fro')^2;
Jpseudo = norm(Anew\Aold-eye(n),'fro')^2;
% Jopt =
%     3.2361
% Jpseudo = 
%     8.0901

% problem 5c
yold = Aold*x;
xhat = B*yold>0.5;
% calculate bit error rate (BER)
BERopt = sum(~x==xhat)/n;
BERpseudo = sum(~x==(Anew\yold>0.5))/n;
% BERopt =
%     0.0333
% BERpseudo =
%     0.1000