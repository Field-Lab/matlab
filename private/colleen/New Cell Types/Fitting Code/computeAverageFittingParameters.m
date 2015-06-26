function [init_mean, init_sd, init_color_weight, init_scale, init_tau, init_n_filters] = computeAverageFittingParameters(fitting_results)
counter = 1;
for i = 1:length(fitting_results)
    if ~isempty(fitting_results{i})
        fit.mean(counter, :) = [fitting_results{i}.center_point_x, fitting_results{i}.center_point_y];
        fit.sd(counter, :) = [fitting_results{i}.center_sd_x, fitting_results{i}.center_sd_y];
        fit.angle(counter) = [fitting_results{i}.center_rotation_angle];
        fit.color_weight(counter, :) = [fitting_results{i}.color_weight_a fitting_results{i}.color_weight_b fitting_results{i}.color_weight_c];
        fit.scale(counter, :) = [fitting_results{i}.scale_one fitting_results{i}.scale_two]; 
        fit.tau(counter,:) = [fitting_results{i}.tau_one fitting_results{i}.tau_two]; 
        fit.n_filters(counter) = [fitting_results{i}.n_filters]; % expand to n1 and n2
        counter = counter+1;
    end
end


init_mean = mean(fit.mean);
init_sd = mean(fit.sd);
init_color_weight = mean(fit.color_weight);
init_scale = mean(fit.scale);
init_tau = mean(fit.tau);
init_n_filters = mean(fit.n_filters);

    