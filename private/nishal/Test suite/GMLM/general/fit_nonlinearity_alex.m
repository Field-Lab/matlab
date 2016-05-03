function [fit_params, gof,fitres1,fitres] = fit_nonlinearity_alex(asr, gs, nbins)

bin_size = ceil(size(gs,2)/nbins);
tmp = sort(gs);
 %gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];
 
% use this if want to fit for saturation!
gs_bins = [tmp(1):(tmp(end)*1.001-tmp(1))/nbins :tmp(end)*1.001];

% gs_bins = [];
% for ibin=0:nbins
% gs_bins = [gs_bins,prctile(gs,100*ibin/nbins)];
% end

nonlinearity=zeros(length(gs_bins)-1,1);
mean_gs_bin=zeros(length(gs_bins)-1,1);
for k=1:length(gs_bins)-1
    nonlinearity(k)=mean(asr(gs>=gs_bins(k) & gs<gs_bins(k+1)));
    mean_gs_bin(k) = mean(gs(gs>=gs_bins(k) & gs<gs_bins(k+1)));
end
mean_gs_bin = mean_gs_bin(~isnan(nonlinearity));
nonlinearity = nonlinearity(~isnan(nonlinearity));
% 

fitres = fit(mean_gs_bin ,nonlinearity, 'a*x./(b+x)','StartPoint',[100, 100], 'Lower', [0.1 0.1 0.1], 'Upper', [2000 2000 2000]);
[fitres1, gof] = fit(gs',asr, 'a*x./(b+x)', 'StartPoint', [fitres.a, fitres.b], 'Lower', [0.1 0.1 0.1], 'Upper', [2000 2000 2000]);
fit_params.a= fitres1.a;
fit_params.b=fitres1.b;

% figure;
plot(mean_gs_bin/120,nonlinearity/120,'b--*','LineWidth',2);
hold on;
plot(mean_gs_bin/120,fitres(mean_gs_bin)/120,'r','LineWidth',2);
hold on;
plot(mean_gs_bin/120,fitres1(mean_gs_bin)/120,'r--','LineWidth',2);

hold on;
histogram(gs(:)/120,'normalization','Probability')


fitres = fit(mean_gs_bin ,nonlinearity, 'scale*normcdf(x,mu,sigma)','StartPoint',[1, max(nonlinearity), 1.5], 'Lower', [0.1 0.1 0.1], 'Upper', [2000 2000 2000]);
[fitres1, gof] = fit(gs',asr, 'scale*normcdf(x,mu,sigma)', 'StartPoint', [fitres.mu, fitres.scale, fitres.sigma], 'Lower', [0.1 0.1 0.1], 'Upper', [2000 2000 2000]);


hold on;
plot(mean_gs_bin/120,fitres(mean_gs_bin)/120,'m','LineWidth',2);
hold on;
plot(mean_gs_bin/120,fitres1(mean_gs_bin)/120,'m--','LineWidth',2);
