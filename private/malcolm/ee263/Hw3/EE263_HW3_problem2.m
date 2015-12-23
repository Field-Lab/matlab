tol = 1e-6;
kmax = 50;
J = zeros(kmax-1,1);
z = ones(N,1);
alpha = zeros(2*m+1,1);
S = ones(N,1);
for i = 1:m
    S = [S (1j*omega).^i]; %#ok<*AGROW>
end
for i = 1:m
    S = [S -h.*(1j*omega).^i];
end

for k = 2:kmax
    htilde = h./z;
    Stilde = diag(1./z) * S;
    Stilde_split = [real(Stilde) ; imag(Stilde)];
    htilde_split = [real(htilde) ; imag(htilde)];
    alpha_prime = Stilde_split \ htilde_split;
    J(k-1) = norm(Stilde * alpha_prime - htilde)/N;
    if norm(alpha_prime-alpha)<tol
        break
    else
        alpha = alpha_prime;
        z = polyval([alpha(17:-1:10);1],1i*omega);
    end
    if k==2
        first_alpha = alpha;
    end
end

% first_alpha =
% 
%     1.2605
%     6.5631
%    36.8476
%    33.6859
%   145.2066
%    70.4987
%   184.2841
%    46.0000
%    74.4660
%    -3.4991
%    30.5001
%   -38.1088
%   176.2871
%   -80.0641
%   319.2128
%   -46.9151
%   177.4030

% alpha =
% 
%    1.0e+03 *
% 
%     0.0014
%     0.0321
%     0.1324
%     0.6485
%     1.7410
%     2.0083
%     4.7318
%     2.0876
%     3.0723
%    -0.0040
%     0.0177
%    -0.1038
%     0.6039
%    -1.0201
%     3.5301
%    -1.2528
%     4.1754

% J(1)
% 
% ans =
% 
%     0.0640
% 
% J(end)
% 
% ans =
% 
%      0

% make plots
subplot(2,1,1);
semilogy(2:kmax,J);
xlabel('k');
ylabel('J');
title('Hw3 P2b');

subplot(2,1,2);
H_first = polyval(first_alpha(9:-1:1),1j*omega)./ ...
    polyval([first_alpha(17:-1:10);1],1j*omega);
H_last = polyval(alpha(9:-1:1),1j*omega)./ ...
    polyval([alpha(17:-1:10);1],1j*omega);
semilogy(omega,abs(H_first),'b',omega,abs(H_last),'r--');
xlabel('omega');
ylabel('abs(H)');