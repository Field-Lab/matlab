function [fit_params, binned_spike_rate, bin_vals] = fit_spike_response(generator, cell_spike_rate, varargin)

p = inputParser;

p.addRequired('generator', @isnumeric);
p.addRequired('cell_spike_rate', @isnumeric);
p.addParamValue('num_pts_per_bin', 200, @isnumeric);
p.addParamValue('verbose', false', @islogical);
p.addParamValue('fig_num', 0, @isnumeric)
p.addParamValue('fit_function', 'parabola', @ischar);

p.parse(generator, cell_spike_rate, varargin{:});


num_pts_per_bin = p.Results.num_pts_per_bin;

[sorted_generator_vals, sorted_generator_indices] = sort(generator,'ascend');

num_bins = floor(length(generator)./num_pts_per_bin);
bin_vals = zeros(num_bins,1);
binned_spike_rate = zeros(num_bins,1);

for bn = 1:num_bins
    begin_bin = 1 + ((bn-1) * num_pts_per_bin);
    end_bin = ((bn-1) * num_pts_per_bin) + num_pts_per_bin;
    temp_gen_indices = sorted_generator_indices(begin_bin:end_bin);
    temp_gen_signals = sorted_generator_vals(begin_bin:end_bin);
    
    mean_of_bin = mean(temp_gen_signals);
    mean_spike_rate_in_bin = mean(cell_spike_rate(temp_gen_indices));
    
    bin_vals(bn) = mean_of_bin;
    binned_spike_rate(bn) = mean_spike_rate_in_bin;
end

switch p.Results.fit_function
    
    case 'parabola'
    
        init_coefs(1) = 1;
        init_coefs(2) = 0.5;
        fit_params = nlinfit(bin_vals, binned_spike_rate, 'parabola_fit', init_coefs);
        fit_values = parabola_fit(fit_params, bin_vals);

    case 'exp'
        
        init_coefs(1) = 1;
        init_coefs(2) = 0.5;
        init_coefs(3) = 0.5;
        fit_params = nlinfit(bin_vals, binned_spike_rate, 'exp_fit', init_coefs);
        fit_values = exp_fit(fit_params, bin_vals);

    otherwise
        
        error('fit_function specification is not recognized')
       
end
        
if p.Results.verbose
    if p.Results.fig_num == 0
        figure
    else
        figure(p.Results.fig_num)
    end
    dot_noise = normrnd(0, 0.05,1, length(cell_spike_rate));
    plot(generator, cell_spike_rate+dot_noise, '.', 'Color', [0.5 0.5 0.5])
    hold on
    plot(bin_vals, binned_spike_rate, 'r.-');
    
    % get fit values and plot
    plot(bin_vals, fit_values, 'b')
    hold off
    
    xlabel('generator signal')
    ylabel('spike count')
    title('gray = data; red = binned data; blue = fit') 
    axis([-1 1 -1 5])
    %legend('data', 'binned data', 'fit')
    %legend boxoff
end

    





    
    