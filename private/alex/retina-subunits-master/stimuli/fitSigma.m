function prsEst = fitSigma(fit1,fit2,testMode)

% set up an fmincon routine
rndMat = randn(fit1.C,fit1.C);
initPrs = [vector(rndMat+rndMat')];

% create handle to objective function
fun = @(prs) fitSigma_errFun(prs(:),fit1,fit2);
x0 = initPrs(:);

if testMode == 1
    DerivCheck_Elts(fun,x0);
    return
end

opts = optimset('Display','iter');
prsEst = fminunc(fun,x0,opts);


%-------------------------------------------
function [err] = fitSigma_errFun(sigma,fit1,fit2)

sigma = reshape(sigma,fit1.C,fit1.C);

%%
term1 = 0;
for is1=1:fit1.S
    for is1_p=1:fit1.S
        term1_1 = fit1.B_s(is1)*fit1.B_s(is1_p);
        vec = (fit1.A_sc(is1,:)+fit1.A_sc(is1_p,:));
        term1_2 = exp(vec*sigma*vec');
        term1_update = term1_1 + term1_2;
        term1 = term1 + term1_update;
    end
end
%%
term2 = 0;
for is2=1:fit2.S
    for is2_p=1:fit2.S
        term2_1 = fit2.B_s(is2)*fit2.B_s(is2_p);
        vec = (fit2.A_sc(is2,:)+fit2.A_sc(is2_p,:));
        term2_2 = exp(vec*sigma*vec');
        term2_update = term2_1 + term2_2;
        term2 = term2 + term2_update;
    end
end
%%
term3 = 0;
for is1=1:fit1.S
    for is2=1:fit2.S
        term3_1 = fit1.B_s(is1)*fit2.B_s(is2);
        vec = (fit1.A_sc(is1,:)+fit2.A_sc(is2,:));
        term3_2 = exp(vec*sigma*vec');
        term3_update = term3_1 + term3_2;
        term3 = term3 + term3_update;
    end
end
%%
err = (term1 + term2 - 2*term3);
