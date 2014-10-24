load fit-2-subunits.mat
fit1 = fit;
%fit1.B_s = ones(size(fit1.B_s))*mean(fit1.B_s);

load fit-8-subunits.mat
fit2 = fit;
%fit2.B_s = ones(size(fit2.B_s))*mean(fit2.B_s);

%%
iCount = 1;
x = zeros(fit_2.C,500);
for i=1:500
    try
        x(:,iCount) = fitStim(fit1,fit2,0);
        iCount = iCount + 1;
    catch
        fprintf('skipping this one');
    end
end
x = x(:,1:iCount-1);

%%
%x = normMat(randn(size(x)),'L2',1);
% eval the output for model 1

    tmp = normMat(randn(size(x(:,i))),'L2',1);
    tmp = U(:,7);
    %tmp = randn(size(x(:,i)));
    %tmp = x(:,i);
%[Y_st_1] = evalNonLin((fit1.I_sc.*fit1.A_sc)*x(:,i),fit1.f);
[Y_st_1] = evalNonLin((fit1.I_sc.*fit1.A_sc)*tmp,fit1.f);
[Z_t_1] = evalNonLin((fit1.J_s.*fit1.B_s)*Y_st_1,fit1.g);

% eval the output for model 2
%[Y_st_2] = evalNonLin((fit2.I_sc.*fit2.A_sc)*x(:,i),fit2.f);
[Y_st_2] = evalNonLin((fit2.I_sc.*fit2.A_sc)*tmp,fit2.f);
[Z_t_2] = evalNonLin((fit2.J_s.*fit2.B_s)*Y_st_2,fit2.g);

foo1 = Z_t_1;
foo2 = Z_t_2;
foo1,
foo2,
%%
for i=1:8
scatter(fit_2.locs_c(:,1),fit_8.locs_c(:,2),50,x(:,i),'filled');
axis off
colormap(cjet);
set(gca,'YDir','reverse');
set(gca,'CLim',[-0.2 0.2]);
drawnow; pause(0.1);
end

%%

%term1 = bsxfun(@times,fit_2.A_sc,fit_2.B_s');
%term2 = bsxfun(@times,fit_8.A_sc,fit_8.B_s');
term1 = fit_2.A_sc;
term1_prime = bsxfun(@times,fit_2.A_sc,fit_2.B_s');
term2 = fit_8.A_sc;
term2_prime = bsxfun(@times,fit_8.A_sc,fit_8.B_s');
%%
load fit-2-subunits.mat
fit1 = fit;

load fit-8-subunits.mat
fit2 = fit;
%%
out = fitSigma(fit1,fit2,0);
%%
[U S V] = svd(cov(x'));
scatter(fit1.locs_c(:,1),fit1.locs_c(:,2),50,U(:,2),'filled');
set(gca,'YDir','reverse');
set(gca,'CLim',[-0.2 0.2]);
axis off; colormap(cjet);

