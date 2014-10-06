function spacing = infer_electrode_spacing(positions)
% INFER_ELECTRODE_SPACING    Use the position coordinates to infer the average distance between electrodes
%
% 2010-06 phli - Abstracted out of plot_ei_.m
%
interpoint_distance_matrix = ipdm(positions, 'subset', 'nearest', 'result', 'struct');
spacing = median(interpoint_distance_matrix.distance);