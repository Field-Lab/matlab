function [stas_sp,newcomputation]=fit_spatial_sta_for_nulling(stas,precomputed,matlab_cellids,CellMask)

% p = inputParser;
% p.addParamValue('center_point_x', 3, @isnumeric);
% p.addParamValue('center_point_y', 3, @isnumeric);
% p.addParamValue('sd_x', 1, @isnumeric);
% p.addParamValue('sd_y', 1, @isnumeric);
% p.addParamValue('amp_scale', 1, @isnumeric);
% p.addParamValue('rotation_angle', 0, @isnumeric);
% p.addParamValue('x_dim', 5, @isnumeric);
% p.addParamValue('y_dim', 5, @isnumeric);
if(~isempty(precomputed.matlab_cellids) & ~(isempty(precomputed.stas_sp)))
    display('Using precomputed');
    I=find(precomputed.matlab_cellids==matlab_cellids);
    stas_sp=precomputed.stas_sp{I(1)};
    newcomputation=0;
    return;
else
    display('Computing fits');
    fit_info=fit_sta_sequence(stas, 'fit_temporal',false,'fit_center',true,'fit_surround',true,'verbose',false);
    stas_sp=make_Gaussian_two_d('center_point_x',fit_info.center_point_x,'center_point_y',fit_info.center_point_y,'sd_x',fit_info.center_sd_x,'sd_y',fit_info.center_sd_y,'amp_scale',fit_info.surround_amp_scale,'rotation_angle',fit_info.center_rotation_angle,'x_dim',fit_info.x_dim,'y_dim',fit_info.y_dim);
    stas_sp=stas_sp.*CellMask; % this step is debatable .. though I feel this is very important!
    newcomputation=1;
end

end