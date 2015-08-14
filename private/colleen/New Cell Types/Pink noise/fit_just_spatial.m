
function [final_fit_params,x ,y] = fit_just_spatial(sta, num_gauss, mark_params)
% 
% sta = datarun.stas.stas{1};
% 
% figure; imagesc(sta(:,:,2,25))


[x,y]=meshgrid(1:size(sta,2),1:size(sta,1));
x=x(:);
y=y(:);

sig_stixels = significant_stixels(sta, 'select','thresh', 'thresh', mark_params.thresh);


params.sig_stixels = sig_stixels;
rf = rf_from_sta(sta, params);
if isempty(rf)
    final_fit_params = nan(6,1);
    return
end

if size(rf,3) > 1
    sta_one = rf(:,:,2);
else
    sta_one = rf;
end

if sum(full(sig_stixels(:))) == 0
    final_fit_params = nan(6,1);
    return
end
          biggestBlob = ExtractNLargestBlobs(full(sig_stixels), 1);
        real_stix = biggestBlob;
[i,j] =find(real_stix);

sta_one = sta_one/max(sta_one(:));


z=sta_one(:);
ellipse_size = num_gauss;

a = std(j);
b = std(i);
h  = mean(j);
k  =mean(i);
angle = 2*pi*rand(1);
A =1;


params = [h;k;a;b;angle;A];

[final_fit_params, fval] = fminsearch(@(params)spatial_fit_error(x,y,z, params),params);

h = final_fit_params(1);
k = final_fit_params(2);
a = ellipse_size*final_fit_params(3);
b = ellipse_size*final_fit_params(4);
angle = final_fit_params(5);
A = final_fit_params(5);

xhat = (x - h)*cos(angle) - (y-k)*sin(angle);
yhat = (x - h)*sin(angle) + (y-k)*cos(angle);
U = (xhat/a).^2 + (yhat/b).^2;
F = A*exp(-U/2);
                    
fit_shaped = reshape(F, size(sta_one,1),size(sta_one,2));
% figure; imagesc(fit_shaped)
%  drawEllipse(h, k, a, b, -angle)
% 
%  axis equal
% 
% figure; imagesc(sta_one)
% 
%  f = drawEllipse(h, k, a, b, -angle);
%  set(f, 'LineWidth', 1.5, 'Color', [0 0 0 ])
%  axis equal


 
 
% 
%  f(x,y) = a1*exp(-(x-x0)^2/(2*sigmax^2)-(y-y0)^2/(2*sigmay^2))
%  
%  [final_fit_params, fval] = fminsearch(@(fit_params)sta_fit_error(sta,...
%                               fit_params, input_params(fixed_indices),...
%                               fit_indices, fixed_indices,...
%                               verbose, mark_params),...
%                               input_params(fit_indices),...
%                               optimset(optim{:}));
%                           
%                           
% % sta_one = sta_one(:,41:80);
% 
% % sig_stixels = significant_stixels(sta_one);
% % [a,b] = find(sig_stixels);
% 
% % y = size(rf,1);
% % x = size(rf,2);
% 
% % x1 = sort(repmat([1:x]', [y, 1]));
% % y1 = repmat([1:y]', [x,1]);
% 
% pairs = [a, b];
% a'*a
% % [zb, ab, bb, alphab] = fitellipse(pairs)
% % A = [x1'*x1 x1'*y1 y1'*x1 y1'*y1]
% % 
% % b = -log(permute(rf(:,:,2), [2,1]))
% % sol = A\b
% % sta_one = squeeze(sta(:,:,2,25));
% 
% % sta_reshape = reshape(sta_one, size(sta_one,1)*size(sta_one,2), size(sta_one,3))
%      % Plot the results
%       plot(pairs(1,:), pairs(2,:), 'ro')
%       hold on
%       plotellipse(zb, ab, bb, alphab, 'b--')
