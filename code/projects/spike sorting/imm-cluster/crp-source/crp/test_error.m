%example to make algorithm crash.
%load included spikes.mat file and run this script.

%initially I thought that increasing the number of datapoints might
%avoid this problem. this turned out not to be the case. it seems as though
%this error is independent of the number of datapoints--it is only
%dependent on the dimensionality of each individual datapoint.

%this is the error message received:
% ??? Error using ==> mvnrnd at 88
% SIGMA must be a symmetric positive semi-definite matrix.
% 
% Error in ==> sampler at 177
%         new_table_mean = mvnrnd(mu_0',new_table_covariance'/k_0)';
% 
% Error in ==> test_sampler at 42
% [class_id, phi, K_record, lP_record, alpha_record] = sampler(spikes', 500);

%should we use five dimensions or three?
FIVE = true;

% running this (each spike is 5 dimensional) will result in an error:
if FIVE
    [class_id, phi, K_record, lP_record, alpha_record] = sampler(spikes', 10000);
end

% running this (each spike is 3 dimensional) will work fine
if ~FIVE
    [class_id, phi, K_record, lP_record, alpha_record] = sampler(spikes(:,1:3)', 500);
end

