
x = data_large(data_rows, 3:end);

y = data_large(data_rows, 2);

% xtest = csvread('data/processed/x_test.csv');
% R = randperm(size(xtest,1),1000);
% xtest = xtest(R,:);
% ytest = csvread('data/processed/y_test.csv');
% 
% ytest = ytest(R,1);
% xtest=[xtest(:,3:11),xtest(:,13)];
% add the intercept term
x = [ones(size(y)), x];
% xtest = [ones(size(ytest)), xtest];
n = size(x, 2);
% numTestDocs = length(ytest);


%% Newton's method
  % run logistic regression using newton's method
  h = @(t, x) 1./(1+exp(- x * t));
  f = @(t) transpose(x) * (y - h(t, x));
  df = @(t) -transpose(x) * diag(h(t, x).*(1-h(t, x))) * x;
  x0 = zeros(n, 1);
  [theta, it] = newton(f, df, x0);
hypo = 1./(1+exp(-x*theta));
hypo(hypo<=0.5) = 0;
hypo(hypo>0.5) = 1;

error=0;
for i=1:numTestDocs
  if (ytest(i) ~= hypo)
    error=error+1;
  end
end

%Print out the classification error on the test set
error = error/numTestDocs
sumHypo = sum(hypo)