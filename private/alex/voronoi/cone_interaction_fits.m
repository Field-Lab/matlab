function [params_ref_all, params_x_all, params_y_all,params_xy_all, ...
    resnorm_x_all, resnorm_y_all, resnorm_xy_all, x_all, y_all] = ...
    cone_interaction_fits(inputs, spikes, center_cones)

spikes(spikes>size(inputs,2)) = [];
spikes_tmp = spikes;
spike_rate=zeros(size(inputs,2),1);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp

%% POLYNOMIAL
acc =[];
for cone1 = 1:length(center_cones)
    cone1
    other_cones = 1:length(center_cones);
    other_cones(cone1)=[];
    for cone2 = other_cones
        
        mean2 = mean(inputs(cone2, :));
        std2 = std(inputs(cone2, :));
        
        inds_mid = find(inputs(cone2, :) > (mean2-std2) & inputs(cone2, :) < (mean2+std2));
        cone_input_mid = inputs(cone1, inds_mid);
        firing_rate_mid = spike_rate(inds_mid);
        
        inds_min = find(inputs(cone2, :) < (mean2-std2));
        cone_input_min = inputs(cone1, inds_min);
        firing_rate_min = spike_rate(inds_min);
        
        inds_max = find(inputs(cone2, :) > (mean2+std2));
        cone_input_max = inputs(cone1, inds_max);
        firing_rate_max = spike_rate(inds_max);
        
        
        figure
        plot(cone_input_min, firing_rate_min, '*')
        hold on
        plot(cone_input_mid, firing_rate_mid, '*')
        plot(cone_input_max, firing_rate_max, '*')
        
        % start points - least square
        [start_points_one, r1, start_points_one2, r2, start_points, r3, start_points_y, r4] = ...
            fit_one_poly(cone_input_min',firing_rate_min, cone_input_max', firing_rate_max);
        
        % log likelihood from x shift fit using least squares
        xscale = start_points(6);
        x_scaled = cone_input_min' + xscale;
        x_com = [x_scaled; cone_input_max'];
        ydata_com = [firing_rate_min; firing_rate_max];
        k = start_points(1)*x_com.^4 + start_points(2)*x_com.^3 + start_points(3)*x_com.^2 + start_points(4)*x_com + start_points(5);
        k(k<=0) = 1e-6;
        y_prelog = k .^ydata_com .*exp(-k) ./ factorial(ydata_com);
        y = -sum(log(y_prelog));
        
        % maximize loglik
        [mleres_one1, mleres_one2, mleres1, mleres2] = ...
            get_mle(cone_input_min',firing_rate_min, cone_input_max', firing_rate_max, ...
            start_points_one, start_points_one2, start_points, start_points_y);
        
        % log likelihood first curve fit maximizing loglik
        x_com = cone_input_min';
        ydata_com = firing_rate_min;
        k = mleres_one1(1)*x_com.^4 + mleres_one1(2)*x_com.^3 + mleres_one1(3)*x_com.^2 + mleres_one1(4)*x_com + mleres_one1(5);
        k(k<=0) = 1e-6;
        y_prelog = k .^ydata_com .*exp(-k) ./ factorial(ydata_com);
        y = -sum(log(y_prelog));
        
        % log likelihood second curve fit maximizing loglik
        x_com = cone_input_max';
        ydata_com = firing_rate_max;
        k = mleres_one2(1)*x_com.^4 + mleres_one2(2)*x_com.^3 + mleres_one2(3)*x_com.^2 + mleres_one2(4)*x_com + mleres_one2(5);
        k(k<=0) = 1e-6;
        y_prelog = k .^ydata_com .*exp(-k) ./ factorial(ydata_com);
        y = -sum(log(y_prelog));
        
        % log likelihood from x shift fit maximizing loglik
        xscale = mleres1(6);
        x_com = [cone_input_min'; cone_input_max' + xscale];
        ydata_com = [firing_rate_min; firing_rate_max];
        k = mleres1(1)*x_com.^4 + mleres1(2)*x_com.^3 + mleres1(3)*x_com.^2 + mleres1(4)*x_com + mleres1(5);
        k(k<=0) = 1e-6;
        y_prelog = k .^ydata_com .*exp(-k) ./ factorial(ydata_com);
        y = -sum(log(y_prelog));
        y_x = y;
        
        % log likelihood from y shift fit maximizing loglik
        x_com = cone_input_min';
        ydata_com = [firing_rate_min; firing_rate_max];
        k1 = mleres2(1)*x_com.^4 + mleres2(2)*x_com.^3 + mleres2(3)*x_com.^2 + mleres2(4)*x_com + mleres2(5);
        x_com = cone_input_max';
        k2 = mleres2(1)*x_com.^4 + mleres2(2)*x_com.^3 + mleres2(3)*x_com.^2 + mleres2(4)*x_com + mleres2(5) + mleres2(6);
        k = [k1; k2];
        k(k<=0) = 1e-6;
        y_prelog = k .^ydata_com .*exp(-k) ./ factorial(ydata_com);
        y = -sum(log(y_prelog));
        
        acc(cone1, cone2)=y_x-y;
        
        % individual fits with x shift from combo fit
        ydata_com = [firing_rate_min; firing_rate_max];
        xscale = mleres1(6);
        x = cone_input_min';
        k1 = mleres_one1(1)*x.^4 + mleres_one1(2)*x.^3 + mleres_one1(3)*x.^2 + mleres_one1(4)*x + mleres_one1(5);
        x1 = cone_input_max' + xscale;
        k2 = mleres_one2(1)*x1.^4 + mleres_one2(2)*x1.^3 + mleres_one2(3)*x1.^2 + mleres_one2(4)*x1 + mleres_one2(5);
        k = [k1; k2];
        k(k<=0) = 1e-6;
        y_prelog = k .^ydata_com .*exp(-k) ./ factorial(ydata_com);
        y = -sum(log(y_prelog));
        
        % individual fits with y shift from combo fit
        ydata_com = [firing_rate_min; firing_rate_max];
        yshift = mleres2(6);
        x = cone_input_min';
        k1 = mleres_one1(1)*x.^4 + mleres_one1(2)*x.^3 + mleres_one1(3)*x.^2 + mleres_one1(4)*x + mleres_one1(5);
        x1 = cone_input_max';
        k2 = mleres_one2(1)*x1.^4 + mleres_one2(2)*x1.^3 + mleres_one2(3)*x1.^2 + mleres_one2(4)*x1 + mleres_one2(5) + yshift;
        k = [k1; k2];
        k(k<=0) = 1e-6;
        y_prelog = k .^ydata_com .*exp(-k) ./ factorial(ydata_com);
        y = -sum(log(y_prelog));
    end
end


figure
plot(acc)
legend(int2str(center_cones))

figure
imagesc(acc)


figure
x=-0.15:0.01:0.15+mleres1(6);
y = mleres1(1)*x.^4 + mleres1(2)*x.^3 + mleres1(3)*x.^2 + mleres1(4)*x + mleres1(5);
plot(x,y, 'linewidth', 2)
hold on
% plot(cone_input_min, firing_rate_min, '+')
% plot(cone_input_max+mleres1(6), firing_rate_max, 'x')

x=-0.15:0.01:0.15;
y1 = mleres_one1(1)*x.^4 + mleres_one1(2)*x.^3 + mleres_one1(3)*x.^2 + mleres_one1(4)*x + mleres_one1(5);
xscale = mleres1(6);
y2 = mleres_one2(1)*x.^4 + mleres_one2(2)*x.^3 + mleres_one2(3)*x.^2 + mleres_one2(4)*x + mleres_one2(5);
plot(x,y1, 'linewidth', 2)
plot(x+xscale,y2, 'linewidth', 2)


figure
x=-0.15:0.01:0.15;
y = mleres2(1)*x.^4 + mleres2(2)*x.^3 + mleres2(3)*x.^2 + mleres2(4)*x + mleres2(5);
plot(x,y, 'linewidth', 2)
hold on
% plot(cone_input_min, firing_rate_min, '+')
% plot(cone_input_max+mleres1(6), firing_rate_max, 'x')

y1 = mleres_one1(1)*x.^4 + mleres_one1(2)*x.^3 + mleres_one1(3)*x.^2 + mleres_one1(4)*x + mleres_one1(5);
yshift = mleres2(6);
y2 = mleres_one2(1)*x.^4 + mleres_one2(2)*x.^3 + mleres_one2(3)*x.^2 + mleres_one2(4)*x + mleres_one2(5) - yshift;
plot(x,y1, 'linewidth', 2)
plot(x,y2, 'linewidth', 2)







figure
x=-0.2:0.01:0.2;
y = mleres2(1)*x.^4 + mleres2(2)*x.^3 + mleres2(3)*x.^2 + mleres2(4)*x + mleres2(5);
plot(x,y)
hold on
plot(cone_input_min, firing_rate_min+mleres2(6), '+')
plot(cone_input_max, firing_rate_max, 'x')


% Poisson
ydata = firing_rate_min';
x = cone_input_min;
y = res.p1*x.^4 + res.p2*x.^3 + res.p3*x.^2 + res.p4*x + res.p5;
y(y<=0) = 1e-6;
y_prelog = y.^ydata .*exp(-y) ./ factorial(ydata);
log(y_prelog);
-sum(log(y_prelog))

mleres = get_mle(cone_input,spike_rate(inds)', res)

ydata = tmp_sr';
x = cone_input;
y = mleres(1)*x.^4 + mleres(2)*x.^3 + mleres(3)*x.^2 + mleres(4)*x + mleres(5);
y(y<=0) = 1e-6;
y_prelog = y.^ydata .*exp(-y) ./ factorial(ydata);
log(y_prelog);
-sum(log(y_prelog))

x=-0.15:0.01:0.15;
y = mleres(1)*x.^4 + mleres(2)*x.^3 + mleres(3)*x.^2 + mleres(4)*x + mleres(5);
plot(x,y)

legend('Bin', 'poly 3', 'poly 4', 'loglik')











cone1 = 2;
cone2 = 5;
mean2 = mean(inputs(cone2, :));
std2 = std(inputs(cone2, :));
% inds = find(inputs(cone2, :) > (mean2-std2) & inputs(cone2, :) < (mean2+std2));
inds = find(inputs(cone2, :) < mean2);
cone_input = inputs(cone1, inds(2:2:end));
tmp_sr = spike_rate(inds);

figure
plot(cone_input, tmp_sr, '*')

% bin first
nbins_cone1 = 20;
tmp = sort(cone_input);
clear cone1_contr
for j=1:nbins_cone1-1
    cone1_contr(j) = tmp(round(size(tmp,2)*j/nbins_cone1));
end
x = [tmp(1); cone1_contr'; tmp(end)];
x = x(1:end-1)+diff(x)/2;
cone1_contr = [-10 cone1_contr 10];
cone1_uncond_rate = zeros(1, nbins_cone1);
cone1_valid = cell(1,nbins_cone1);
cone1_std = 0;
for j = 1:nbins_cone1
    cone1_valid{j} = find(cone_input>cone1_contr(j) & cone_input<=cone1_contr(j+1));
    cone1_uncond_rate(j) = mean(tmp_sr(cone1_valid{j}));
    cone1_std(j) = std(tmp_sr(cone1_valid{j}))/sqrt(length(cone1_valid{j}));
end
figure
plot(x, cone1_uncond_rate)
hold on
plot(x, cone1_uncond_rate+2*cone1_std)
plot(x, cone1_uncond_rate-2*cone1_std)
figure
plot(x, cone1_std)

x=-0.15:0.01:0.15;
[res gof]=  fit(cone_input',tmp_sr, 'poly3');
y = res.p1*x.^3 + res.p2*x.^2 + res.p3*x + res.p4;
plot(x,y)

[res gof] =  fit(cone_input',tmp_sr, 'poly4');
y = res.p1*x.^4 + res.p2*x.^3 + res.p3*x.^2 + res.p4*x + res.p5;
plot(x,y)

% orig
% x = cone_input;
% y = res.p1*x.^4 + res.p2*x.^3 + res.p3*x.^2 + res.p4*x + res.p5;
% y(y<=0) = 1e-6;
% negloglik = -sum(log(1 - y))

% Poisson
ydata = tmp_sr';
x = cone_input;
y = res.p1*x.^4 + res.p2*x.^3 + res.p3*x.^2 + res.p4*x + res.p5;
y(y<=0) = 1e-6;
y_prelog = y.^ydata .*exp(-y) ./ factorial(ydata);
log(y_prelog);
-sum(log(y_prelog))

mleres = get_mle(cone_input,spike_rate(inds)', res)

ydata = tmp_sr';
x = cone_input;
y = mleres(1)*x.^4 + mleres(2)*x.^3 + mleres(3)*x.^2 + mleres(4)*x + mleres(5);
y(y<=0) = 1e-6;
y_prelog = y.^ydata .*exp(-y) ./ factorial(ydata);
log(y_prelog);
-sum(log(y_prelog))

x=-0.15:0.01:0.15;
y = mleres(1)*x.^4 + mleres(2)*x.^3 + mleres(3)*x.^2 + mleres(4)*x + mleres(5);
plot(x,y)

legend('Bin', 'poly 3', 'poly 4', 'loglik')



% 
% log_y = (ydata.*log(y) - y) - log(factorial(ydata));
% -sum(log_y)




custpdf = @(x, p1, p2, p3, p4, p5)(p1*x.^4 + p2*x.^3 + p3*x.^2 + p4*x + p5);

get_mle(x,res)
mle(tmp_sr, 'pdf', custpdf, 'start', [res.p1 res.p2 res.p3 res.p4 res.p5])


x = cone_input;
predic = res.p1*x.^4 + res.p2*x.^3 + res.p3*x.^2 + res.p4*x + res.p5;
etol = 10^-6;
tmp_fr =tmp_sr;
tmp_fr(tmp_fr==0) = etol;
predic(predic<=etol) = etol;

% compute liklihood
loglik = -(sum(tmp_fr(:).*log(predic(:))) - sum(predic(:)));

% lower half
inds = find(inputs(cone2, :)<mean(inputs(cone2, :)));
cone_input = inputs(cone1, inds);
tmp_sr = spike_rate(inds)/max(spike_rate(inds));

figure
hold on

[res gof]=  fit(cone_input',tmp_sr, 'poly3');
y = res.p1*x.^3 + res.p2*x.^2 + res.p3*x + res.p4;
plot(x,y)

[res gof] =  fit(cone_input',tmp_sr, 'poly4');
y = res.p1*x.^4 + res.p2*x.^3 + res.p3*x.^2 + res.p4*x + res.p5;
plot(x,y)
x = cone_input;
predic = res.p1*x.^4 + res.p2*x.^3 + res.p3*x.^2 + res.p4*x + res.p5;
etol = 10^-6;
tmp_fr =tmp_sr;
tmp_fr(tmp_fr==0) = etol;
predic(predic<=etol) = etol;

loglik = -(sum(tmp_fr(:).*log(predic(:))) - sum(predic(:)));

% joint
inds = find(inputs(cone2, :)>mean(inputs(cone2, :)));
cone_input = inputs(cone1, inds);
tmp_sr = spike_rate(inds)/max(spike_rate(inds));

inds1 = find(inputs(cone2, :)<mean(inputs(cone2, :)));
cone_input1 = inputs(cone1, inds);
tmp_sr1 = spike_rate(inds)/max(spike_rate(inds));

predic = res.p1*x.^4 + res.p2*x.^3 + res.p3*x.^2 + res.p4*x + res.p5;

[res gof] =  fit(cone_input',tmp_sr, 'poly4');



%% OLD BINNING



nbins_cone1 = 6;

for cone1 = 1:length(center_cones)
    cone1
    other_cones = 1:length(center_cones);
    other_cones(cone1) = [];
    
    
    tmp = sort(inputs(cone1,:));
    
    clear cone1_contr
    for j=1:nbins_cone1-1
        cone1_contr(j) = tmp(round(size(tmp,2)*j/nbins_cone1));
    end
    x = [tmp(1); cone1_contr'; tmp(end)];
    x = x(1:end-1)+diff(x)/2;
    
    cone1_contr = [-10 cone1_contr 10];
    cone1_uncond_rate = zeros(1, nbins_cone1);
    cone1_valid = cell(1,nbins_cone1);
    cone1_std = 0;
    for j = 1:nbins_cone1
        cone1_valid{j} = find(inputs(cone1,:)>cone1_contr(j) & inputs(cone1,:)<=cone1_contr(j+1));
        cone1_uncond_rate(j) = mean(spike_rate(cone1_valid{j}));
        cone1_std(j) = std(spike_rate(cone1_valid{j}));
    end
    
    
    datx = x;
    
    for cone2 = other_cones
        y2 = cone1_uncond_rate';
        tmp = sort(inputs(cone2,:));
        cone2_valid = find(inputs(cone2,:)>tmp(ceil(size(inputs,2)/5*4)));
        cond_rate = zeros(nbins_cone1,1);
        for j = 1:nbins_cone1
            tt = intersect(cone1_valid{j}, cone2_valid);
            cond_rate(j) = mean(spike_rate( tt(1:2:end)));
            cond_std1(j) = std(spike_rate(tt(1:2:end)));
        end
        y2 = [y2 cond_rate];
        
        cone2_valid = find(inputs(cone2,:)<tmp(ceil(size(inputs,2)/5*2)));
        cond_rate = zeros(nbins_cone1,1);
        for j = 1:nbins_cone1
            tt = intersect(cone1_valid{j}, cone2_valid);
            cond_rate(j) = mean(spike_rate( tt(1:2:end)));
            cond_std2(j) = std(spike_rate(tt(1:2:end)));
        end
        y2 = [y2 cond_rate];
        y2(isinf(y2)) = 0;
        
        x2 = repmat(datx, 1,size(y2,2));
        
%         figure;plot(x2,y2)
        
        eps2 = [cone1_std; cond_std1; cond_std2]';

        [params_ref, params_x, params_y, params_xy, resnorm_x, resnorm_y, resnorm_xy] = fit_2_curves(x2,y2, eps2);
        
        
        params_ref_all{cone1,cone2} = params_ref;
        params_x_all{cone1,cone2} = params_x;
        params_y_all{cone1,cone2} = params_y;
        params_xy_all{cone1,cone2} = params_xy;
        resnorm_x_all(cone1, cone2) = resnorm_x;
        resnorm_y_all(cone1, cone2) = resnorm_y;
        resnorm_xy_all(cone1, cone2) = resnorm_xy;
        x_all{cone1} = x2;
        y_all{cone1,cone2} = y2;
    end
    
end


% 
% 
% figure
% nbins_cone1 = 50;
% for cone1 = 1:length(center_cones)
%     cone1
%     other_cones = 1:length(center_cones);
%     other_cones(cone1) = [];
%     
%     
%     tmp = sort(inputs(cone1,:));
%     
%     clear cone1_contr
%     for j=1:nbins_cone1-1
%         cone1_contr(j) = tmp(round(size(tmp,2)*j/nbins_cone1));
%     end
%     x = [tmp(1); cone1_contr'; tmp(end)];
%     x = x(1:end-1)+diff(x)/2;
%     
%     cone1_contr = [-10 cone1_contr 10];
%     cone1_uncond_rate = zeros(1, nbins_cone1);
%     cone1_valid = cell(1,nbins_cone1);
%     cone1_std = 0;
%     for j = 1:nbins_cone1
%         cone1_valid{j} = find(inputs(cone1,:)>cone1_contr(j) & inputs(cone1,:)<=cone1_contr(j+1));
%         cone1_uncond_rate(j) = mean(spike_rate(cone1_valid{j}));
%         cone1_std(j) = std(spike_rate(cone1_valid{j}));
%     end
%     
%     plot(x, cone1_uncond_rate)
%     tm = (max(cone1_uncond_rate)- min(cone1_uncond_rate))/4;
% 
%     
%     % min value
%     clear new_contr
%     thresh = mean(cone1_uncond_rate)-tm;
%     tt = find(cone1_uncond_rate<thresh, 1, 'last')
%     new_contr(1) = cone1_contr(tt);
%     thresh = mean(cone1_uncond_rate);
%     tt = find(cone1_uncond_rate<thresh, 1, 'last')
%     new_contr(2) = cone1_contr(tt);
%     thresh = mean(cone1_uncond_rate)+tm;
%     tt = find(cone1_uncond_rate<thresh, 1, 'last')
%     new_contr(3) = cone1_contr(tt);
%     
%     x = [tmp(1); new_contr'; tmp(end)];
%     x = x(1:end-1)+diff(x)/2;
%     
%     nbins_cone1 = 4
%     new_contr = [-10 new_contr 10];
%     cone1_uncond_rate = zeros(1, nbins_cone1);
%     cone1_valid = cell(1,nbins_cone1);
%     cone1_std = 0;
%     for j = 1:nbins_cone1
%         cone1_valid{j} = find(inputs(cone1,:)>new_contr(j) & inputs(cone1,:)<=new_contr(j+1));
%         cone1_uncond_rate(j) = mean(spike_rate(cone1_valid{j}));
%         cone1_std(j) = std(spike_rate(cone1_valid{j}));
%     end
%     
%     plot(x,cone1_uncond_rate)
% 
%     
%     
%     cond_rate = zeros(nbins_cone1,1);
%     for j = 1:nbins_cone1
%         tt = intersect(cone1_valid{j}, cone2_valid);
%         cond_rate(j) = mean(spike_rate( tt(1:2:end)));
%         cond_std1(j) = std(spike_rate(tt(1:2:end)));
%     end
%     y2 = [y2 cond_rate];
%     
%     
%     datx = x;
%     
% %     for cone2 = other_cones
%     for cone2 = 5
%         y2 = cone1_uncond_rate';
%         tmp = sort(inputs(cone2,:));
%         
%         clear cone2_contr
%         nbins_cone2 = 16;
%         for j=1:nbins_cone2-1
%             cone2_contr(j) = tmp(round(size(tmp,2)*j/nbins_cone2));
%         end
%         x = [tmp(1); cone2_contr'; tmp(end)];
%         x = x(1:end-1)+diff(x)/2;
%         
%         cone2_contr = [-10 cone2_contr 10];
%         cone2_uncond_rate = zeros(1, nbins_cone2);
%         cone2_valid = cell(1,nbins_cone2);
%         for j = 1:nbins_cone2
%             cone2_valid{j} = find(inputs(cone2,:)>cone2_contr(j) & inputs(cone2,:)<=cone2_contr(j+1));
%             cone2_uncond_rate(j) = mean(spike_rate(cone2_valid{j}));
%         end
%         clear cone2_valid
%         tm = (max(cone2_uncond_rate)- min(cone2_uncond_rate))/4;
%         
%         thresh = mean(cone2_uncond_rate)-tm;
%         fr = 0; i = 1;
%         while  (fr< thresh  | i<7) & i<26
%             cone2_valid = find(inputs(cone2,:)<tmp(1000*i));
%             fr = mean(spike_rate(cone2_valid));
%             i = i+1;
%         end        
%         cond_rate = zeros(nbins_cone1,1);
%         for j = 1:nbins_cone1
%             tt = intersect(cone1_valid{j}, cone2_valid);
%             cond_rate(j) = mean(spike_rate( tt(1:2:end)));
%             cond_std1(j) = std(spike_rate(tt(1:2:end)));
%         end
%         y2 = [y2 cond_rate];
%         
%         thresh = mean(cone2_uncond_rate)+tm;
%         fr = 1; i = 1;
%         while (fr> thresh  | i<7) & i<26
%             cone2_valid = find(inputs(cone2,:)>tmp(end-1000*i));
%             fr = mean(spike_rate(cone2_valid));
%             i = i+1;
%         end        
%         cond_rate = zeros(nbins_cone1,1);
%         for j = 1:nbins_cone1
%             tt = intersect(cone1_valid{j}, cone2_valid);
%             cond_rate(j) = mean(spike_rate( tt(1:2:end)));
%             cond_std1(j) = std(spike_rate(tt(1:2:end)));
%         end
%         y2 = [y2 cond_rate];
%         
%         
% %         cone2_valid = find(inputs(cone2,:)<tmp(ceil(size(inputs,2)/5*2)));
% %         cond_rate = zeros(nbins_cone1,1);
% %         for j = 1:nbins_cone1
% %             tt = intersect(cone1_valid{j}, cone2_valid);
% %             cond_rate(j) = mean(spike_rate( tt(1:2:end)));
% %             cond_std2(j) = std(spike_rate(tt(1:2:end)));
% %         end
% %         y2 = [y2 cond_rate];
% %         y2(isinf(y2)) = 0;
%         
%         x2 = repmat(datx, 1,size(y2,2));
%         
%         subplot(3,3,cone1);
%         hold on
%         plot(x2,y2)
%         
%     end
%     
% end
% %
% %
% % 
% % figure;
% plot(cone1_uncond_rate)
% hold on
% figure
% for nbins_cone2 = 4:2:16
%     tmp = sort(inputs(cone2,:));
%     clear cone2_contr
%     for j=1:nbins_cone2-1
%         cone2_contr(j) = tmp(round(size(tmp,2)*j/nbins_cone2));
%     end
%     x = [tmp(1); cone2_contr'; tmp(end)];
%     x = x(1:end-1)+diff(x)/2;
%     
%     cone2_contr = [-10 cone2_contr 10];
%     cone2_uncond_rate = zeros(1, nbins_cone2);
%     cone2_valid = cell(1,nbins_cone2);
%     for j = 1:nbins_cone2
%         cone2_valid{j} = find(inputs(cone1,:)>cone2_contr(j) & inputs(cone1,:)<=cone2_contr(j+1));
%         cone2_uncond_rate(j) = mean(spike_rate(cone2_valid{j}));
%     end
%     hold on
%     plot(x,cone2_uncond_rate)
% end
% line([-0.1 0.1], [mean(spike_rate), mean(spike_rate)], 'color', 'k')
% size(cone2_valid{1},2)
% 
% 
% tm = (max(cone2_uncond_rate)-min(cone2_uncond_rate))/10;
% 
% thresh = mean(spike_rate)+tm*2.5;
% tmp = sort(inputs(cone2,:), 'descend');
% td = 1;
% i = 1;
% while td>thresh
%     contr_val = tmp(2000*i);
%     cone2_valid = find(inputs(cone2,:)>=contr_val);
%     td = mean(spike_rate(cone2_valid));
%     i = i+1;
% end
% 
% cond_rate = zeros(nbins_cone1,1);
% for j = 1:nbins_cone1
%     cond_rate(j) = mean(spike_rate(intersect(cone1_valid{j}, cone2_valid)));
%     cond_std1(j) = std(spike_rate(intersect(cone1_valid{j}, cone2_valid)));
% end
% plot(cond_rate)
%         
% y2 = [cone1_uncond_rate' cond_rate];
% eps2 = [cone1_std; cond_std1];
% 
% thresh = mean(spike_rate)-tm*1.5;
% tmp = sort(inputs(cone2,:), 'ascend');
% td = -1;
% i = 1;
% while td<thresh
%     contr_val = tmp(2000*i);
%     cone2_valid = find(inputs(cone2,:)<=contr_val);
%     td = mean(spike_rate(cone2_valid));
%     i = i+1;
% end
% 
% cond_rate = zeros(nbins_cone1,1);
% for j = 1:nbins_cone1
%     cond_rate(j) = mean(spike_rate(intersect(cone1_valid{j}, cone2_valid)));
%     cond_std1(j) = std(spike_rate(intersect(cone1_valid{j}, cone2_valid)));
% end
% plot(cond_rate)
% 
% y2 = [y2 cond_rate];
% eps2 = [eps2; cond_std1]';
%      
% [params_ref, params_x, params_y, params_xy, resnorm_x, resnorm_y, resnorm_xy] = fit_2_curves(x2,y2, eps2);
%        
% 
% 
% % %plot
% % x = datx(:,1);
% % sat   = g_t(1);
% % sigma = g_t(2);
% % mu = g_t(3);
% % sh = g_t(4);
% % y = sat .* normcdf(x, mu, sigma)+sh;
% % figure
% % hold on
% % plot(x,y, '-*')
% % plot(x,raw_dat(:,1))
% 
% 
% % % plot
% % x = [-1:0.05:0.5];
% % n = size(datx,2);
% % sat   = p(n+1);
% % sigma = p(n+2);
% % mu = 0;
% % sh = p(n+3);
% % y = sat .* normcdf(x, mu, sigma)+sh;
% % figure
% % hold on
% % plot(x,y, '-*')
% % hold on
% % plot(datx(:,1)+p(1),raw_dat(:,1))
% % plot(datx(:,1)+p(2),raw_dat(:,2))
% % plot(datx(:,1)+p(3),raw_dat(:,3))
% 
% 
% % % plot
% % x = datx(:,1);
% % n = size(datx,2);
% % sat   = g(n+1);
% % sigma = g(n+2);
% % mu = g(n+3);
% % y = sat .* normcdf(x, mu, sigma);
% % figure
% % hold on
% % plot(x,y, '-*')
% % plot(x,raw_dat(:,1)-g(1))
% % plot(x,raw_dat(:,2)-g(2))
% % plot(x,raw_dat(:,3)-g(3))