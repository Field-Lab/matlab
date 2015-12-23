% problem 4a
% define target vector etilde
etilde = [reshape(eye(m,m),m^2,1); zeros(n*m,1)];
% calculate xtilde for each value of lambda
lambda = 4:2:16;
xtilde = nan(m*n,length(lambda));
J1 = nan(length(lambda),1);
J2 = nan(length(lambda),1);
for i = 1:length(lambda)
    disp(i);
    PhiTilde = makePhiTilde(Phi,m,n,lambda(i));
    xtilde(:,i) = PhiTilde\etilde;
    J = (PhiTilde*xtilde(:,i)-etilde).^2;
    J1(i) = sum(J(1:m^2));
    J2(i) = sum(J(m^2+1:end))/lambda(i);
end

% optimal J when lambda=10
% >> J1(lambda==10)+10*J2(lambda==10)
% ans =
%     8.3068

% plot tradeoff curve, J2 vs J1, and mark point where lambda=10
plot(J2,J1,'-.k');
hold on;
plot(J2(lambda==10),J1(lambda==10),'x','MarkerSize',20);
title('Problem 4a-1');
xlabel('J2');
ylabel('J1');

% make stem plot of w(x(15))
figure();
stem(1:m,Phi*xtilde(14*n+1:15*n,lambda==10),'k.');
grid on;
title('Problem 4a-2');
xlabel('i');
ylabel('w(x^{(15)})_i');

% problem 4b
% use x corresponding to lambda=10
x = reshape(xtilde(1:m*n,lambda==10),n,m);
% define matrix Wtilde
Wtilde = nan(k*m,m);
for i = 1:m
    Wtilde((i-1)*k+1:i*k,:) = repmat((Phi*x(:,i))',k,1);
end
% define target vector Ytilde
Ytilde = reshape(Y',m*k,1);
zhat = Wtilde\Ytilde;

% problem 4c
zhat1 = (Phi*x)'\Y(:,1);
Erms = sqrt(sum(ztrue-zhat).^2/m);
Erms1 = sqrt(sum(ztrue-zhat1).^2/m);
% Erms =
%     0.0972
% Erms1 =
%     0.7126 

% plot ztrue, zhat, and zhat1
plot(1:50,ztrue,'-k',1:50,zhat,'--xk',1:50,zhat1,':ok');
grid off;
title('Problems 4b and c');
xlabel('i');
ylabel('z_i');
legend({'ztrue','zhat','zhat1'},'location','northwest')