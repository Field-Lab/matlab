function SNR = get_cone_weight_SNR(mosaic_weights, selection)


num_rgcs = size(mosaic_weights,2);

SNR = zeros(num_rgcs,1);
for cc = 1:num_rgcs
    temp_weights = mosaic_weights(:,cc);
    
    % normalize weights
    temp_weights = temp_weights ./ max(temp_weights);
    
    temp_rb_sd = robust_std(temp_weights);
    
%    temp_cone_indices = find(selection(:,cc));
    sig_cone_weights = temp_weights(selection(:,cc));
    mean_sig_cone_weights = mean(sig_cone_weights);
    
    temp_SNR = mean_sig_cone_weights ./ temp_rb_sd;
    
    SNR(cc) = temp_SNR;
    
end