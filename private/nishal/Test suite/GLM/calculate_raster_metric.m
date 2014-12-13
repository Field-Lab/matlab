function response = calculate_raster_metric(response,metric_params)

if(strcmp(metric_params.type,'psth-sd'))
    % Careful about effects of bin size!
    response.psth_sd = sqrt(var(mean(response.spksGen,2)));
end

end