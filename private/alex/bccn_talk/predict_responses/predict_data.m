function predicted_rate = predict_data(sta, gensig_bins, nonlinearity, test_inputs)

% filter raw inputs with STA
filt_inputs=zeros(size(test_inputs,1),size(test_inputs,2)-size(sta,2)+1);
for i=1:size(test_inputs,1)
    filt_inputs(i,:)=conv(test_inputs(i,:),sta(i,:),'valid');
end
gen_sig=sum(filt_inputs,1);
gen_sig=gen_sig/max(gen_sig);

% pass generator signal through nonlinearity
predicted_rate = zeros(size(gen_sig));
for j=1:size(gensig_bins,2)-1
    a = gen_sig >= gensig_bins(j) & gen_sig < gensig_bins(j+1);
    predicted_rate(a)= nonlinearity(j);
end
