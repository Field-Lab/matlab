function ttf_out = transform_tf(ttf,scale)
% scale
%ttf_out = (ttf-mean(ttf))*scale + mean(ttf);

% scale
% r = 20;
% ttf_upsample = interp(ttf,r);
% iidx=1:r:r*length(ttf);
% iidx = ceil(iidx/scale);
% ttf_out = ttf_upsample(iidx);

% delay ? 
r = 20;
ttf_upsample = interp(ttf,r);
iidx=1:r:r*length(ttf);
iidx = max(ceil(iidx)-scale,1);
ttf_out = ttf_upsample(iidx);

figure;
plot(ttf);hold on;plot(ttf_out);legend('Original','Stretched');
pause(0.1);
end