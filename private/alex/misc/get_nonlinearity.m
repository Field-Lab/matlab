function [nonlinearity, gs_bins] = get_nonlinearity(data, gs, nbins)

% get bins    
bin_size = ceil(size(gs,2)/nbins);
tmp = sort(gs);
gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];

nonlinearity=zeros(length(gs_bins)-1,1);
for k=1:length(gs_bins)-1
    nonlinearity(k)=mean(data(gs>=gs_bins(k) & gs<gs_bins(k+1)));
end