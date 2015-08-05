init_space = reshape(spatialfilterfromSTA(fittedGLM.cellinfo.WN_STA, fittedGLM.rawfit.ROIcoord.xvals, fittedGLM.rawfit.ROIcoord.yvals), [13,13]);
init_SU = 0.1*ones(3);
init_SU(5)=0.4;

figure;

subplot(2,5,1)
imagesc(init_SU)
colormap gray
caxis([0 0.5])
axis image
axis off
title('Initial')
for iter = 1:4
subplot(2,5,iter+1)
imagesc(reshape(fittedGLM.rawfit.iter{iter}.SU, [3 3]))
colormap gray
caxis([0 0.5])
axis image
axis off
title(['Iter' num2str(iter)])
end

subplot(2,5,6)
imagesc(init_space)
colormap gray
caxis([-0.8 1])
axis image
axis off
for iter = 1:4
subplot(2,5,iter+6)
imagesc(reshape(fittedGLM.rawfit.iter{iter}.nonSU(fittedGLM.rawfit.paramind.space1), [13 13]))
colormap gray
caxis([-0.8 1])
axis image
axis off
end
