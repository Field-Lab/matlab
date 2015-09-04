function predicted_rate = convert_gs(gs, gs_bins, nonlinearity)

predicted_rate = zeros(size(gs));
for j=1:size(gs_bins,2)-1
    a = gs >= gs_bins(j) & gs < gs_bins(j+1);
    predicted_rate(a)= nonlinearity(j);
end