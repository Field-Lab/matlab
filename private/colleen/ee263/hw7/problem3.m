% problem 3

% create matrix Ctilde (this method is a little janky but works)
D = blkdiag(C,C,C,C,C,C,C,C,C,C);
Ctilde = blkdiag(D,D,D,D,D,D,D,D,D,D);

% initialize Btilde
Btilde = B*v(1);
for i = 2:100
    Bcurr = B*v(i);
    for j = 1:i-1
        Bcurr = Bcurr + (A^j)*B*v(i-j);
    end
    Btilde = [Btilde; Bcurr];
end

% calculate singular values of Ctilde*Btilde
[~,s,q] = svd(Ctilde*Btilde);
qmax = q(:,1);
qmin = q(:,3);

% calculate associated values of D
Dmax = 0.1*norm(Ctilde*Btilde*qmax);
Dmin = 0.1*norm(Ctilde*Btilde*qmin);