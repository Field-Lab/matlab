% problem 2
close all; clc;

% 2a
% create unit ball
theta = 0:0.1:2*pi;
z = [cos(theta); sin(theta)];
figure();
% calculate both ellipse representations for all 4 measurements
B = {[];[];[];[]};
xhat = {[];[];[];[]};
S = {[];[];[];[]};
for i = 1:M
    xhat{i} = A{i}\y{i};
    yhat = A{i}*xhat{i};
    c = sqrt(alpha(i)^2-norm(y{i}-yhat)^2);
    [~,sigma,V] = svd(A{i},0);
    B{i} = c*V*diag(1./diag(sigma));
    S{i} = (1/c^2)*A{i}'*A{i};
    ellipse = B{i}*z+repmat(xhat{i},1,length(theta));
    if i==1
        plot(ellipse(1,:),ellipse(2,:));
    elseif i==2
        plot(ellipse(1,:),ellipse(2,:),'x');
    elseif i==3
        plot(ellipse(1,:),ellipse(2,:),'o');
    else
        plot(ellipse(1,:),ellipse(2,:),'--');
    end
    hold on;
end
title('Problem 2a');

% 2b
areas = nan(M,1);
for i = 1:M
    areas(i) = pi*prod(svd(B{i}));
end
% areas =
%    53.1580
%    39.1085
%    55.2141
%    49.1166

% 2c
areas_hat = nan(M,1);
for i = 1:M
    sigma = svd(B{i});
    % generate square with side 2*sigma1 centered at xhat
    w1 = 2*max(sigma)*(rand(N,1)-0.5);
    w2 = 2*max(sigma)*(rand(N,1)-0.5);
    w = [w1 w2]+repmat(xhat{i}',N,1);
    % check membership of points
    phat = 0;
    for j = 1:N
        phat = phat + ((w(j,:)-xhat{i}')*S{i}*(w(j,:)-xhat{i}')' <= 1);
    end
    phat = phat/N;
    areas_hat(i) = phat*(2*max(sigma))^2;
end
% areas_hat =
%    52.6573
%    39.8521
%    55.4253
%    48.8122

RMSE = sqrt(sum((areas-areas_hat).^2)/M);
% RMSE =
%     0.4850

intersection_areas = nan(M,M);
for i = 1:M
    for j = i+1:M
        % generate large rectangle
        phat = 0;
        sigma_i = max(svd(B{i}));
        sigma_j = max(svd(B{j}));
        xhati = xhat{i}; xhatj = xhat{j};
        xhat_mean = (xhati + xhatj)/2;
        d1 = max([abs(xhati(1)+sigma_i) ; abs(xhati(1)-sigma_i) ;...
                  abs(xhatj(1)+sigma_j) ; abs(xhatj(1)-sigma_i)]);
        d2 = max([abs(xhati(2)+sigma_i) ; abs(xhati(2)-sigma_i) ;...
                  abs(xhatj(2)+sigma_j) ; abs(xhatj(2)-sigma_i)]);
        w1 = 2*d1*(rand(N,1)-0.5);
        w2 = 2*d2*(rand(N,1)-0.5);
        w = [w1 w2]+repmat(xhat_mean',N,1);
        % check membership of points
        phat = 0;
        for k = 1:N
            phat = phat + (...
                ((w(k,:)-xhat{i}')*S{i}*(w(k,:)-xhat{i}')' <= 1) & ...
                ((w(k,:)-xhat{j}')*S{j}*(w(k,:)-xhat{j}')' <= 1) );
        end
        phat = phat/N;
        intersection_areas(i,j) = phat*4*d1*d2;
    end
end
% intersection_areas =
%        NaN   13.8043   10.7333   14.4438
%        NaN       NaN   13.0810   37.1368
%        NaN       NaN       NaN   16.3619
%        NaN       NaN       NaN       NaN