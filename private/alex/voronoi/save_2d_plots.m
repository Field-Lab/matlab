function save_2d_plots(xdata1, xdata2, y, mllparams_x, mllparams_y, cone1, cone2, cellID, LLR, matr_LLR)


figure;
set(gcf, 'position', [-1886         164        1048         930]);

clear fr
wei = mean([mllparams_x(end), mllparams_y(end)]);
xdata1_tmp = xdata1*wei;
% tt = min(max(xdata1_tmp), max(xdata2));
% tt1 = max(min(xdata1_tmp), min(xdata2));
% tmp1 = linspace(tt1, tt, 7);
% tmp2 = tmp1;
tmp1 = linspace(min(xdata1_tmp), max(xdata1_tmp), 10);
tmp2 = linspace(min(xdata2), max(xdata2), 10);
for i = 2:length(tmp1)
    for j = 2:length(tmp2)
        inds = find(xdata1_tmp>tmp1(i-1) & xdata1_tmp<=tmp1(i) & xdata2>tmp2(j-1) & xdata2<=tmp2(j));
        fr(i-1,j-1) = mean(y(inds));
    end
end
[x1_new, x2_new] = meshgrid(tmp1(1:end-1),tmp2(1:end-1));

subplot(2,2,3)
surface(x1_new, x2_new, fr');
view(3)
xlabel(['cone ',int2str(cone1) ,' weight ', num2str(wei)])
ylabel(['cone ',int2str(cone2)])
zlabel('firing rate')
title(['LLR ', num2str(LLR)], 'Interpreter', 'none')
path2save = ['/Volumes/Analysis/2016-02-17-4/subunits/surface_plots/cell_', int2str(cellID),'/'];


all_array = zeros(7,7);
all_array(2:6,2:6) = reshape(matr_LLR,5,5);
all_array = all_array/max(abs(all_array(:)))/2+0.5;
all_array = repmat(all_array,1,1,3);
tmp = all_array;
x0 = 0.5;
k = 10;
tmp = ones(size(tmp))./(1+exp(-k*(tmp-x0)));
subplot(2,2,4)
imagesc(tmp);
m = mean(matr_LLR(:));
s = std(matr_LLR(:));

if m/s>0
    title({['cone ',int2str(cone1) ,' and ',int2str(cone2)], ...
        [num2str(m), ' / ', num2str(s), '  =  ', num2str(m/s)], ...
        'separate subunits likelier'});
else
        title({['cone ',int2str(cone1) ,' and ',int2str(cone2)], ...
        [num2str(m), ' / ', num2str(s), '  =  ', num2str(m/s)], ...
        'same subunit likelier'});
end
axis off
set(gca,'dataaspectratio', [1 1 1])


if 1
    p = mllparams_x;
    
%     tt = min(max(xdata1), max(xdata2));
%     tt1 = max(min(xdata1), min(xdata2));
%     tmp = linspace(tt1, tt, 15);
    
    x1 = x1_new;
    x2 = x2_new;
    
    sat   = p(1);
    k = p(2);
    mu = p(3);
    c = p(4);
    a = p(5);
    y_interim = ones(size(x1)).*sat./(1+ exp(-k*(a*x1 + x2 - mu)))+c;
    
%     sat   = p(1);
%     sigma = p(2);
%     mu = p(3);
%     sh = p(4);
%     a = p(5);
%     y_interim = sat .* normcdf(a*x1 + x2, mu, sigma) + sh;
    subplot(2,2,1)
    surface(x1, x2, y_interim);
    view(3)
    title('x shift')

    
    p = mllparams_y;
    
    sat   = p(1);
    k = p(2);
    mu = p(3);
    c = p(4);
    a = p(5);
    y_interim = ones(size(x1)).*sat./(1+ exp(-k*(a*x1 - mu))) + ones(size(x1)).*sat./(1+ exp(-k*(x2 - mu))) + c;
%     sat   = p(1);
%     sigma = p(2);
%     mu = p(3);
%     sh = p(4);
%     a = p(5);
%     y_interim = sat .* (normcdf(a*x1, mu, sigma) + normcdf(x2, mu, sigma)) + sh;
    
    subplot(2,2,2)
    surface(x1, x2, y_interim);
    view(3)
    title('y shift')
end

% if ~isdir(path2save)
%     mkdir(path2save);
% end
saveas(gcf, [path2save, 'cone_',int2str(cone1) ,'_',int2str(cone2), '.bmp']);
close all

% saveas(gcf, ['/Volumes/Analysis/2011-12-13-2/cone_data/manual/surface_plots_3736/data_only/data_', ...
%     int2str(cone1), '_',int2str(cone2),'.svg' ]);
% close all
