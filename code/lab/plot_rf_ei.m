function plot_rf_ei(datarun,cell_id)
% plot_rf_ei     plot an RF and overlay the EI
%
% usage:  plot_rf_ei(datarun,cell_id)
%
% arguments:     datarun - datarun struct
%                cell_id - cell id
%
%
% 2010-04  gauthier
%





% BODY OF THE FUNCTION

% plot RF
plot_rf(datarun,cell_id,'array',true)

% plot EI
hold on;
plot_ei(datarun,cell_id,'coordinates','sta','pretty_axes',0,'alpha',0,...
    'neg_color','b','pos_color','r','scale',1)

