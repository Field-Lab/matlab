%% Script for validating the gradient using the findif utility

basepars.nvar = length(p0);

pars.basepars = basepars;
pars.stimpars = stimpars;
pars.trainpars = trainpars;
pars.fgname = 'test_hess';

vartest = [1 5 100 150 160 1+144+40+4 1+144+40+60 1+2*44+40+4 1+2*144+40+32  length(p0)-5]

randp0 = 0.1.*randn(size(p0));
for j=vartest
    pars.hesstest_j = j;
    findif(randp0,pars);
    pause(0.25);
    1;
end




%% Test the Hessian function
%p = 0.01.*randn(size(p0));
%p0 = xstar;
basepars.nvar = length(p0);
error2 = zeros(basepars.nvar);

func = @(x) ll_func2(x,basepars,stimpars,trainpars);

[f ll_grad ll_hess] = func(p0);
%ll_hess =  ll_hess_mult(Hinfo,eye(length(p0)),basepars,stimpars,trainpars);

vartest = 1:10;

for i=vartest
    
  % Perturb the ith variable and observe effect the gradient
        
  if (mod(i,1) == 0)
      fprintf('i=%d\n',i);
  end
  res = 16;
  v = zeros(res,basepars.nvar);
  for j=res:res
    p2 = p0; p2(i) = p2(i) + 2^(-j);
    
    [lg_p2 ll_grad2] = func(p2);
     v(j,:) = ll_grad2'; %C(:,j) = cifs2;
  end
  grad_approx = (v - repmat(ll_grad',res,1))./repmat((2.^(-1:-1:-res))',1,basepars.nvar); % res x nvars
  error2(i,:) = ll_hess(i,:) - grad_approx(res,:);
%  if (error > 10^(-3))
%     fprintf('Error of %f in component %d\n',error,i);        
%  end
end