
iframe = 1;
coneFrame = train.X_ct(:,iframe);
scatter(dat.locs_c(:,1),dat.locs_c(:,2),50,coneFrame,'filled');

mnValX = min(dat.locs_c(:,1));
mxValX = max(dat.locs_c(:,1));
mnValY = min(dat.locs_c(:,2));
mxValY = max(dat.locs_c(:,2));
res = 31;
xRange = linspace(mnValX-1,mxValX+1,res);
yRange = linspace(mnValY-1,mxValY+1,res);
coneVar = 0.001;
[X Y] = meshgrid(xRange,yRange);
[distance diff] = getDistances(dat.locs_c,[X(:) Y(:)]);
diff = reshape(permute(diff,[1 3 2]),[dat.nConesFit*res^2,2]);
cones = evalGaussian(diff,[1 0; 0 1]*coneVar);
cones = reshape(cones,[dat.nConesFit res^2]);
cones = normMat(cones,'L2',2);
M = cones;

foo = M'*train.X_ct;
%%
colormap(gray);
for i=1:25
    imagesc(xRange,yRange,reshape(foo(:,i),res,res),[-4 4]);
    %hold on
    %scatter(dat.locs_c(:,1),dat.locs_c(:,2),100,train.X_ct(:,i),'filled');
    axis equal;
    pause(0.1);
    drawnow;
end
%%
ft = zeros(res^2,size(foo,2));
for i=1:size(foo,2)
    stimFFT = fft2(reshape(foo(:,i),res,res));
    stimFFTmag = abs(stimFFT(:));
    ft(:,i) = stimFFTmag;
end
%%
A = reshape(ft.^2,[res^2 size(foo,2)]);
%b = regress(train.R_t',[ones(1,size(foo,2)); A]');
b = ridge((train.R_t)',A',100,0);
%%
%tmpSpk = ft(:,logical(train.R_t));
allSpk = ft.^2;

tmpSpk = bsxfun(@times,allSpk,(train.R_t));

estFiltFFT = sqrt(fftshift(reshape(mean(tmpSpk,2),res,res)));
estFilt = real(fftshift(ifft2(ifftshift(estFiltFFT))));
%%
estFiltFFT = (fftshift(reshape(b(2:end),res,res)));
estFilt = real(fftshift(ifft2(ifftshift(estFiltFFT))));
