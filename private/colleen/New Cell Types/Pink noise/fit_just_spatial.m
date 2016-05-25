
function [final_fit_params,x ,y] = fit_just_spatial(sta, num_gauss, mark_params)

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



