function fit_error = MDLP_error(sta_fit_params, cell_spike_rate, kernel_matrix, data_range)


% the parameters of the nonlinear associated with the STA -- use a half
% squaring nonlinearity with x and y offsets
sta_fit_fun = @(sta_fit_params)sta_output_error(sta_fit_params, kernel_matrix(1,:), cell_spike_rate);

[fit_params, f_val] = fminsearch(sta_fit_fun, sta_fit_params);

% get the LNP model ouput firing rate


num_dims = size(kernel_matrix,1);


if isempty(data_range);
    clear data_range
    data_range(1) = 1;
    data_range(2) = length(cell_spike_rate);
end    
    
fit_spike_rate = cell_spike_rate(data_range(1):data_range(2));

NL_output = zeros(1,length(fit_spike_rate));
for dm = 1:num_dims
    % get NL output from the STA
    if dm == 1
        sta_nl_output = kernel_matrix(1,:);
        sta_nl_output(sta_nl_output < 0) = sta_nl_y_offset;
        sta_nl_output(sta_nl_output >= 0)=  (weights(dm) * (kernel_matrix(1,:)-sta_nl_x_offset).^2) + sta_nl_y_offset;
    else
        NL_output = NL_output + (weights(dm) .* kernel_params(dm).kernel_generator(data_range(1):data_range(2)).^2);
    end
end
%NL_output(NL_output<0) = 0;
fit_error = sqrt(mean((cell_spike_rate - NL_output).^2));



