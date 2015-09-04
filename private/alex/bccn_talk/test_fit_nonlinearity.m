function [nonlinearity, gs_bins] = test_fit_nonlinearity(asr, gs, nbins)

frac = 0.5;
res = zeros(1,20);
for i = 1:20
    a = randperm(length(gs));
    
    % fitting
    inds = a(1:ceil(length(gs)*frac));
    
    gs1 = gs(inds);
    asr1 = asr(inds);
    bin_size = ceil(size(gs1,2)/nbins);
    tmp = sort(gs1);
    gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];
    nonlinearity=zeros(length(gs_bins)-1,1);
    for k=1:length(gs_bins)-1
        nonlinearity(k)=mean(asr1(gs1>=gs_bins(k) & gs1<gs_bins(k+1)));
    end
    [fitres, gof] = fit(gs_bins(1:end-1)',nonlinearity, 'a*normcdf(x,mu,sigma)','StartPoint',[0.1, 1, 1.5]);
    [fitres1, gof1] = fit(gs1',asr1, 'a*normcdf(x,mu,sigma)', 'StartPoint', [fitres.a, fitres.mu, fitres.sigma]);
    % testing
    inds = a(ceil(length(gs)*frac)+1:end);
    gs2 = gs(inds);
    asr2 = asr(inds);
    tt = fitres1.a*normcdf(gs2,fitres1.mu, fitres1.sigma);
%     figure
%     plot(tt)
%     hold on
%     plot(asr2)
    res(i) = corr(tt', asr2);
end
mean(res)
figure
plot(res)



for fr=0:4
    
    a = round([(length(gs)*fr/5+1):length(gs) 1:length(gs)*fr/5-1]);
    % fitting
    inds = a(1:ceil(length(gs)*frac));
    gs1 = gs(inds);
    asr1 = asr(inds);
    bin_size = ceil(size(gs1,2)/nbins);
    tmp = sort(gs1);
    gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];
    nonlinearity=zeros(length(gs_bins)-1,1);
    for k=1:length(gs_bins)-1
        nonlinearity(k)=mean(asr1(gs1>=gs_bins(k) & gs1<gs_bins(k+1)));
    end
    [fitres, gof] = fit(gs_bins(1:end-1)',nonlinearity, 'a*normcdf(x,mu,sigma)','StartPoint',[max(nonlinearity), 1, 1.5]);
    [fitres1, gof1] = fit(gs1',asr1, 'a*normcdf(x,mu,sigma)', 'StartPoint', [fitres.a, fitres.mu, fitres.sigma]);
    % testing
    inds = a(ceil(length(gs)*frac)+1:end);
    gs2 = gs(inds);
    asr2 = asr(inds);
    tt = fitres1.a*normcdf(gs2,fitres1.mu, fitres1.sigma);
    gg(fr+1) = corr(tt', asr2)
end

mean(gg)

gs1 = gs;
asr1 = asr;
bin_size = ceil(size(gs1,2)/nbins);
tmp = sort(gs1);
gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];
nonlinearity=zeros(length(gs_bins)-1,1);
for k=1:length(gs_bins)-1
    nonlinearity(k)=mean(asr1(gs1>=gs_bins(k) & gs1<gs_bins(k+1)));
end
[fitres, gof] = fit(gs_bins(1:end-1)',nonlinearity, 'a*normcdf(x,mu,sigma)','StartPoint',[max(nonlinearity), 1, 1.5]);
[fitres1, gof1] = fit(gs1',asr1, 'a*normcdf(x,mu,sigma)', 'StartPoint', [fitres.a, fitres.mu, fitres.sigma]);
% testing
tt = fitres1.a*normcdf(gs1,fitres1.mu, fitres1.sigma);
corr(tt', asr)



gs1 = gs;
asr1 = asr;
bin_size = ceil(size(gs1,2)/nbins);
tmp = sort(gs1);
gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];
nonlinearity=zeros(length(gs_bins)-1,1);
for k=1:length(gs_bins)-1
    nonlinearity(k)=mean(asr1(gs1>=gs_bins(k) & gs1<gs_bins(k+1)));
end
predicted_rate = zeros(size(gs));
for j=1:size(gs_bins,2)-1
    a = gs >= gs_bins(j) & gs < gs_bins(j+1);
    predicted_rate(a)= nonlinearity(j);
end
corr(predicted_rate', asr)


frac = 0.5
for fr=0:4
    
    a = round([(length(gs)*fr/5+1):length(gs) 1:length(gs)*fr/5-1]);
    % fitting
    inds = a(1:ceil(length(gs)*frac));
    gs1 = gs(inds);
    asr1 = asr(inds);
    bin_size = ceil(size(gs1,2)/nbins);
    tmp = sort(gs1);
    gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];
    nonlinearity=zeros(length(gs_bins)-1,1);
    for k=1:length(gs_bins)-1
        nonlinearity(k)=mean(asr1(gs1>=gs_bins(k) & gs1<gs_bins(k+1)));
    end

    % testing
    inds = a(ceil(length(gs)*frac)+1:end);
    gs2 = gs(inds);
    asr2 = asr(inds);
    predicted_rate = zeros(size(gs2));
    for j=1:size(gs_bins,2)-1
        a = gs2 >= gs_bins(j) & gs2 < gs_bins(j+1);
        predicted_rate(a)= nonlinearity(j);
    end
    gg(fr+1) = corr(predicted_rate', asr2)
end
mean(gg)



frac = 0.1;
res = zeros(1,20);
for i = 1:20
    a = randperm(length(gs));
    
   % fitting
    inds = a(1:ceil(length(gs)*frac));
    gs1 = gs(inds);
    asr1 = asr(inds);
    bin_size = ceil(size(gs1,2)/nbins);
    tmp = sort(gs1);
    gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];
    nonlinearity=zeros(length(gs_bins)-1,1);
    for k=1:length(gs_bins)-1
        nonlinearity(k)=mean(asr1(gs1>=gs_bins(k) & gs1<gs_bins(k+1)));
    end

    % testing
    inds = a(ceil(length(gs)*frac)+1:end);
    gs2 = gs(inds);
    asr2 = asr(inds);
    predicted_rate = zeros(size(gs2));
    for j=1:size(gs_bins,2)-1
        a = gs2 >= gs_bins(j) & gs2 < gs_bins(j+1);
        predicted_rate(a)= nonlinearity(j);
    end
    res(i) =  corr(predicted_rate', asr2);
end

mean(res)
figure
plot(res)
