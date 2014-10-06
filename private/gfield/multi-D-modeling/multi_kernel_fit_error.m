function fit_error = multi_kernel_fit_error(weights, offsets, generator_matrix, spike_rate)

output_sigs = multi_kernel_output(offsets, generator_matrix);

predicted_spike_rate = weights * output_sigs;

fit_error = sqrt(mean((predicted_spike_rate - spike_rate).^2));

% figure(1); clf; hold on
% plot(generator_matrix(1,:), output_sigs(1,:), '.r')
% plot(generator_matrix(1,:), spike_rate, '.k')
% drawnow

