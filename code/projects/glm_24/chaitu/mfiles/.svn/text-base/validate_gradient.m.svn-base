%% Script for validating the gradient using the findif utility

pars2.basepars = basepars;
pars2.stimpars = stimpars;
pars2.trainpars = trainpars;
pars2.fgname = 'llfunc_findif_wrapper';
findif(0.1.*randn(size(p0)),pars2);

%% Test the gradient function
% Analyze the each component of the gradient

nvar = length(p0);


%func2 = @(x) func(x,[]);
func2 =(@(p)ll_func2(p,basepars,stimpars,trainpars)); % @(p)ll_func2_hessmult_precon(p,basepars,stimpars,trainpars,Sinvtotal);%@(x) ll_func2(x,basepars,stimpars,trainpars,[]);
[lg_p ll_grad] = func2(p0);

error = zeros(nvar,1);
% 1 26 66 91 131 156 196 221 261 271
vartest = [1:10];


for i1=vartest
    i1
    if (mod(i1,10) == 0)
        fprintf('i=%d\n',i1);
    end
    res = 16;
    v = zeros(res,1);
    for j1=1:5:res
        p2 = p0; p2(i1) = p2(i1) + 2^(-j1);
        
        [lg_p2] = func2(p2);
        
        %[lg_p2 cifs2 kx2] = train_ll3(p2,basepars,stimpars,trainpars);
        v(j1) = lg_p2; %C(:,j) = cifs2;
    end
    grad_approx = (v - lg_p)'./(2.^(-1:-1:-res));
    error(i1) = abs(ll_grad(i1) - grad_approx(floor(res)));
    error(i1)
    if (error(i1) > 10^(-3))
        fprintf('Error of %f in component %d\n',error(i1),i1);        
    end
end

