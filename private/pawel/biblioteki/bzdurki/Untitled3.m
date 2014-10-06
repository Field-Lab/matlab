fo = fitoptions('method','NonlinearLeastSquares','Normalize','off');
%ok_ = isfinite(AllMovies) & isfinite(Effic);
%if ~all( ok_ )
%    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
%        'Ignoring NaNs and Infs in data' );
%end
%st_ = [0.28946898665270315 0.5608943345498274 ];
%set(fo_,'Startpoint',st_);

ft = fittype('A*(exp(-(x-tau).^2/(2*sigma^2)))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'A', 'tau', 'sigma'});
ft = fittype(gaussian)
t=[0:1:600]
% Fit this model using new data
set(fo,'Startpoint',[400 90 5]);
cf = fit(t',p',ft,fo);
g=cf(t);
figure(4)
plot(t,p,'b-',t,g,'r-')

figure(3)
subplot(1,2,1)
plot(t,p)
subplot(1,2,2)
plot(t,g)