% problem 1b
n = 2;
dmax = 3;
c = nan(2*dmax+1,dmax); % matrix to store answers
for d = 1:dmax
    alpha = eye(2*d+1);
    alpha = alpha(:,n+1); % desired output (1 if m=n, else 0)
    A = nan(2*d+1, 2*d+1);
    dseq = -d:d;
    for m = 0:2*d
        A(m+1,:) = (dseq.^m/factorial(m))';
    end
    c(1:2*d+1,d) = A\alpha; % invert A to compute coefficients c
end

% c =
% 
%     1.0000   -0.0833    0.0111
%    -2.0000    1.3333   -0.1500
%     1.0000   -2.5000    1.5000
%        NaN    1.3333   -2.7222
%        NaN   -0.0833    1.5000
%        NaN       NaN   -0.1500
%        NaN       NaN    0.0111

% problem 1c
delta = 10.^(-6:-1);
f = @(x) exp(cos(10*x));
% actual 2nd deriv
f2 = @(x) 100*(sin(10*x).^2-cos(10*x)).*exp(cos(10*x)); 

 % estimate 2nd deriv using coefficients from part b
f2hat = nan(3,length(delta));
for d = 1:3
    dseq = -d:d;
    for i = 1:length(delta)
        f2hat(d,i) = (1/delta(i)^n)*sum(c(1:2*d+1,d).*f(1+dseq*delta(i))');
    end
end

% make plot
x = 1./delta;
y = abs(f2hat-f2(1));
loglog(x,y(1,:),'-o',x,y(2,:),'--x',x,y(3,:),':^')
xlim([1 1e7])
xlabel('-log(\delta)')
ylabel('log(|f^{(2)}_{hat}(1)-f^{(2)}(1)|)')
legend({'d=1','d=2','d=3'})