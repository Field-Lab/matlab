%spatial_fit_error

function fit_error = spatial_fit_error(x, y, z, params)
h = params(1);

k = params(2);
a = params(3);
b = params(4);
angle = params(5);
A = params(6);

xhat = (x - h)*cos(angle) - (y-k)*sin(angle);
yhat = (x - h)*sin(angle) + (y-k)*cos(angle);
U = (xhat/a).^2 + (yhat/b).^2;
F = A*exp(-U/2);

fit_error = sqrt(mean(( z  - F).^2));


end


% function to plot spatial fit ellispe
function plot_spatial_sd(params)

    % center
    ctr = [params(1) params(2)];
    rad_sd = [params(4) params(3)];
    rot_angle = params(5);
    
    % get points of an ellipse with these parameters
    [X, Y] = drawEllipse([ctr rad_sd (-1*rot_angle)]);
    plot(X, Y, 'k')
    
    % surround
    if params(13) ~= 0
        rad_sd_sur = rad_sd * params(13);
        [X, Y] = drawEllipse([ctr rad_sd_sur (-1*rot_angle)]);
        plot(X, Y, 'r')
    end
end



% function to constrain parameters
function all_params = apply_constraints(all_params)

    % ensure surround_radius >= center_radius
    if all_params(13) < 1; all_params(13) = 1; end
    
    % ensure surround_scale >= 0
    if all_params(14) < 0; all_params(14) = 0; end
    
    % ensure surround_scale < 1
    if all_params(14) >= 1; all_params(14) = 1; end

    % ensure center_radius >= 0
    if all_params(3) < 0; all_params(3) = realmin; end

    % ensure center_radius >= 0
    if all_params(4) < 0; all_params(4) = realmin; end
    
    % ensure the color triplet has unit length (norm = 1)
%     norm_factor = norm([all_params(6), all_params(7), all_params(8)]);
%     all_params(6) = all_params(6) ./ norm_factor;
%     all_params(7) = all_params(7) ./ norm_factor;
%     all_params(8) = all_params(8) ./ norm_factor;
%    
end

