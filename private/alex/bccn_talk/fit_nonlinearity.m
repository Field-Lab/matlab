function [fit_params, gof] = fit_nonlinearity(asr, gs, nbins)

bin_size = ceil(size(gs,2)/nbins);
tmp = sort(gs);
gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];
nonlinearity=zeros(length(gs_bins)-1,1);
for k=1:length(gs_bins)-1
    nonlinearity(k)=mean(asr(gs>=gs_bins(k) & gs<gs_bins(k+1)));
end
fitres = fit(gs_bins(1:end-1)',nonlinearity, 'scale*normcdf(x,mu,sigma)','StartPoint',[1, max(nonlinearity), 1.5], 'Lower', [0.1 0.1 0.1], 'Upper', [15 15 15]);
[fitres1, gof] = fit(gs',asr, 'scale*normcdf(x,mu,sigma)', 'StartPoint', [fitres.mu, fitres.scale, fitres.sigma], 'Lower', [0.1 0.1 0.1], 'Upper', [15 15 15]);

fit_params.scale = fitres1.scale;
fit_params.mu = fitres1.mu;
fit_params.sigma = fitres1.sigma;

% 
% figure
% plot(gs_bins(1:end-1)',nonlinearity)
% hold on
% plot(gs_bins(1:end-1)', fitres.scale*normcdf(gs_bins(1:end-1),fitres.mu,fitres.sigma))

