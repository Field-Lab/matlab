G = nan(N-M+1, M+1);
G(:,1) = my_gamma(u(M:N));
for i = 1:M
    G(:,i+1) = z(M-i+1:N-i+1);
end
y = my_mu(u(M:N));
alpha = G\y;
mse = norm(G*alpha-y)^2/(N-M+1);
 
% alpha =
% 
%    -0.5004
%     0.9949
%     0.1980
%    -0.3513
%    -0.3400
%    -0.0624
%     0.1655
%     0.1503
%     0.0107
%    -0.0745
%    -0.0503
% 
% mse =
% 
%     0.0012

% make plots
subplot(2,1,1)
stem(alpha(2:end))
title('HW3 P3b')

% now calculate u_hat
F = nan(N, M);
for i = 1:M
    F(:,i) = [zeros(i-1,1) ; z(1:N-i+1)];
end
v_hat = F*alpha(2:end);
u_hat = my_psi(v_hat,alpha(1));

% make plots
t = 1:N;
subplot(2,1,2)
plot(t,u,'-',t,z,'--',t,u_hat,'--.',t,u-u_hat,'o')