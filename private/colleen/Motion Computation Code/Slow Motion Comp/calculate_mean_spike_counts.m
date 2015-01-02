% Calculate mean spike counts for each datarun
    %
    % cell_indices1 = get_cell_indices(datarun{2}, {'On parasol'});
    % numspikes1 = zeros(length(cell_indices1),1);
    % for i=1:length(cell_indices1)
    %     numspikes1(i) = length(datarun{2}.spikes{cell_indices1(i)});
    % end
    % 
    % cell_indices2 = get_cell_indices(datarun{2}, {'On midget'});
    % numspikes2 = zeros(length(cell_indices2),1);
    % for i=1:length(cell_indices2)
    %     numspikes2(i) = length(datarun{2}.spikes{cell_indices2(i)});
    % end
    % 
    % if strcmp(run_opt.data_set, '2007-03-27-1')
    %     mean_spikes_onp_exp1(run_opt.data_run-13) = mean(numspikes1);
    %     mean_spikes_onm_exp1(run_opt.data_run-13) = mean(numspikes2);
    % elseif strcmp(run_opt.data_set, '2007-08-24-4')
    %     mean_spikes_onp_exp2(run_opt.data_run-3) = mean(numspikes1);
    %     mean_spikes_onm_exp2(run_opt.data_run-3) = mean(numspikes2);
    % end