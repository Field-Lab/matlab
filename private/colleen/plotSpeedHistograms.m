%Plot histograms for speed estimates for ON parasol cells run 18 all  four
%configs
clear
addpath('/Users/vision/Desktop/GitHub code repository/private/colleen/colleenResults/2007-03-27-1')
estimate(:,1) = importdata('On parasol_data_run_18_config_1.mat')
estimate(:,2) = importdata('On parasol_data_run_18_config_2.mat')
estimate(:,3) = importdata('On parasol_data_run_18_config_3.mat')
estimate(:,4) = importdata('On parasol_data_run_18_config_4.mat')

%1: dark bar, x_delta= 8
%2 dark bar, x_delta = -8
%3 light bar, x_delta = 8
%4 light bar, x_delta = -8

subplot(2,2,1)
histfit(estimate(:,1),20)
title('Dark Bar moving right')

subplot(2,2,2)
histfit(estimate(:,2),20)
title('Dark Bar moving left')

subplot(2,2,3)
histfit(estimate(:,3),20)
title('Bright bar moving right')

subplot(2,2,4)
histfit(estimate(:,4),20)
title('Bright bar moving left')