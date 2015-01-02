function [cell_types, velocity_lim, config_num, trial_estimate_start, tol] = auto_set_params(data_set, data_run)

% Function to autoset parameters for motion scripts
% Malcolm Campbell 2014
% malcolmc@stanford.edu

    if strcmp(data_set, '2007-03-27-1')
        cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};
        if data_run == 12
            velocity_lim = 150;
            config_num = 1;
            trial_estimate_start = 110;
            tol = 1e-3;
        elseif data_run == 13
            velocity_lim = 300;
            config_num = 1;
            trial_estimate_start = 203;
            tol = 1e-4;
        elseif data_run == 14
            velocity_lim = 50;
            config_num = 1;
            trial_estimate_start = 14.6;
            tol = 1e-2;
        elseif data_run == 15
            velocity_lim = 150;
            config_num = 3;
            trial_estimate_start = 110;
            tol = 1e-3;
        elseif data_run == 16
            velocity_lim = 300;
            config_num = 1;
            trial_estimate_start = 203;
            tol = 1e-4;
        elseif data_run == 17
            velocity_lim = 50;
            config_num = 2;
            trial_estimate_start = 14.6;
            tol = 1e-2;
        elseif data_run == 18
            velocity_lim = 150;
            config_num = 1;
            trial_estimate_start = 110;
            tol = 1e-3;
        elseif data_run == 19
            velocity_lim = 300;
            config_num = 1;
            trial_estimate_start = 203;
            tol = 1e-4;
        end
    elseif strcmp(data_set, '2007-08-24-4')
        cell_types = {'Off parasol', 'On midget', 'On parasol'};
        if data_run == 2
            velocity_lim = 150;
            config_num = 2;
            trial_estimate_start = 100;
            tol = 1e-3;
        elseif data_run == 3
            velocity_lim = 300;
            config_num = 3;
            trial_estimate_start = 200;
            tol = 1e-4;
        elseif data_run == 4
            velocity_lim = 300;
            config_num = 1;
            trial_estimate_start = 200;
            tol = 1e-4;
        elseif data_run == 5
            velocity_lim = 300;
            config_num = 2;
            trial_estimate_start = 200;
            tol = 1e-4;
        elseif data_run == 6
            velocity_lim = 150;
            config_num = 3;
            trial_estimate_start = 100;
            tol = 1e-3;
        elseif data_run == 7
            velocity_lim = 150;
            config_num = 1;
            trial_estimate_start = 100;
            tol = 1e-3;
        elseif data_run == 8
            velocity_lim = 50;
            config_num = 1;
            trial_estimate_start = 12.5;
            tol = 1e-2;
        elseif data_run == 9
            velocity_lim = 50;
            config_num = 1;
            trial_estimate_start = 12.5;
            tol = 1e-2;
        elseif data_run == 10
            velocity_lim = 150;
            config_num = 2;
            trial_estimate_start = 50;
            tol = 1e-2;
        elseif data_run == 11
            velocity_lim = 150;
            config_num = 1;
            trial_estimate_start = 50;
            tol = 1e-2;
        end
    end

end