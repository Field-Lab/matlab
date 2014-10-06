function plot_rf_outlines_(sta_fit,color);
%modified for data from PLOT_RF_OUTLINES shlens 2006-07-25
if nargin<3
    color='k';
end

% radius of STA fits in standard standard deviations
radius = 1;


% circle samples
circle_samples = 0:0.1:2*pi;
x_circle = cos(circle_samples);
y_circle = sin(circle_samples);


% grab parameters
sd    = sta_fit.sd;
mn    = sta_fit.mean;
angle = sta_fit.angle;

% rotate by angle and stretch
R = rmatrix2(angle / (2*pi) * 360);
L = radius * [sd(1) 0; 0 sd(2)];
X = R * L * [x_circle; y_circle];
X(:,end+1) = X(:,1);


p = plot(X(1,:)+mn(1), -X(2,:)+mn(2), color);
