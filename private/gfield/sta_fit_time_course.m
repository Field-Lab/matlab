function fit_time_course = sta_fit_time_course(fit_params, fit_type)






switch fit_type
    
    case 'matlab'
        scale_one = fit_params.scale_one;
        scale_two = fit_params.scale_two;
        tau_one = fit_params.tau_one;
        tau_two = fit_params.tau_two;
        n_filters = fit_params.n_filters;

        color_triplet = [fit_params.color_weight_a,...
                fit_params.color_weight_b,...
                fit_params.color_weight_c];

        frame_number = fit_params.frame_number;
        num_colors = fit_params.num_colors;
                    
    case 'obvius'
        scale_one = fit_params.scale(1);
        scale_two = fit_params.scale(2);
        tau_one = fit_params.tau(1);
        tau_two = fit_params.tau(2);
        n_filters = fit_params.n(1);
        
        
        color_triplet = fit_params.color;
        
        frame_number = 29;
        num_colors = length(unique(color_triplet));
        
    otherwise
        error('fit_type must be obvius or matlab')
        
end
            
t_points = 0:1:(frame_number-1);


t_filter_one = scale_one .* (t_points ./ tau_one).^n_filters .* exp(-n_filters*((t_points ./ tau_one) - 1));
t_filter_two = scale_two .* (t_points ./ tau_two).^n_filters .* exp(-n_filters*((t_points ./ tau_two) - 1));
tc = t_filter_one + t_filter_two;
tc = tc(frame_number:-1:1);

fit_time_course = zeros(frame_number, num_colors);
for clr = 1:num_colors
    fit_time_course(:,clr) = tc * color_triplet(clr);
end

