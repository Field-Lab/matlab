[histx histy resSUB resSING coneids weights] = getUDfiles('2011-10-25-5-data001-0',2583);

normOpt = 4;
[resp predic] = getUD(histx,histy,coneids,resSUB.fit_SUB,1);
subplot(2,4,1);
imagesc(resp,[0 1]); axis image off;
subplot(2,4,2);
imagesc(predic,[0 1]); axis image off;
title(sprintf('R2: %.02g',getR2(predic(:),resp(:))));
[resp predic] = getUD(histx,histy,coneids,resSUB.fit_LN,1);
subplot(2,4,4);
imagesc(predic,[0 1]); axis image off;
title(sprintf('R2: %.02g',getR2(predic(:),resp(:))));
[resp predic] = getUD(histx,histy,coneids,resSUB.fit_SUB,normOpt);
subplot(2,4,5);
imagesc(resp,[0 1]); axis image off;
subplot(2,4,6);
imagesc(predic,[0 1]); axis image off;
title(sprintf('R2: %.02g',getR2(predic(:),resp(:))));
[resp predic] = getUD(histx,histy,coneids,resSUB.fit_LN,normOpt);
subplot(2,4,8);
imagesc(predic,[0 1]); axis image off;
title(sprintf('R2: %.02g',getR2(predic(:),resp(:))));

%load('/Users/freemanj11/Dropbox/Projects/RetinaSubunits/anal/2012-04-13-1/data002/subunit/3002-subunit-singletons-0.33-fit.mat')

[resp predic] = getUD(histx,histy,coneids,resSING.fit_SUB,1);
subplot(2,4,3);
imagesc(predic,[0 1]); axis image off;
title(sprintf('R2: %.02g',getR2(predic(:),resp(:))));
[resp predic] = getUD(histx,histy,coneids,resSING.fit_SUB,normOpt);
foo1 = [foo1,resp(:)];
foo2 = [foo2,predic(:)];
subplot(2,4,7);
imagesc(predic,[0 1]); axis image off;
title(sprintf('R2: %.02g',getR2(predic(:),resp(:))));

%coneids = [1 4 7 3]; % for 2011-10-25-5-data003
