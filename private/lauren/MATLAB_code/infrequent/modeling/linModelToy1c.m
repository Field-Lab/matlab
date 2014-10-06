%linearity toy model
clear all

xPlotLim = [0 1.2];
yPlotLim = [0 1.2];

%% constructing the toy model

lam_c1 = -0.45; %cell one, secondary electrode 1
lam_c2 = -0.2; %cell one, secondary electrode 2
lam_c3 = -0.05; %primary electrode different from c1, c2

tp_c1 = 0.6;
tp_c2 = 0.55;
tp_c3 = 0.8;
w1 = 5;
w2 = 2.5;
w3 = 4;

alpha = -3:0.1:3; %secondary:primary ratios

blue = [27 117 187]/255;
red = [190 30 45]/255;
green = [41 180 115]/255;

%% primary electrode + secondary electrode
figure('position', [100 100 200 200])
hold on
t_new_c1 = tp_c1*exp(alpha*lam_c1); %cell 1
t_new_c2 = tp_c2*exp(alpha*lam_c2); %cell 2
t_new_c3 = tp_c3*exp(alpha*lam_c3);

% plot(alpha.*t_new_c1, t_new_c1, '--', 'LineWidth', 2, 'Color', red)
% plot(alpha.*t_new_c2, t_new_c2, '--', 'LineWidth', 2, 'Color', blue)
% plot(t_new_c3, alpha.*t_new_c3, '--', 'LineWidth', 2, 'Color', green) %axes switched because primary/secondary electrodes swapped
% 
% %plot standard deviations
% plot(alpha.*t_new_c1*(1-1/w1), t_new_c1*(1-1/w1), '--', 'Color', red)
% plot(alpha.*t_new_c1*(1+1/w1), t_new_c1*(1+1/w1), '--', 'Color', red)
% 
% plot(alpha.*t_new_c2*(1-1/w2), t_new_c2*(1-1/w2), '--', 'Color', blue)
% plot(alpha.*t_new_c2*(1+1/w2), t_new_c2*(1+1/w2), '--', 'Color', blue)
% 
% plot(t_new_c3*(1-1/w3), alpha.*t_new_c3*(1-1/w3), '--', 'Color', green)
% plot(t_new_c3*(1+1/w3), alpha.*t_new_c3*(1+1/w3), '--', 'Color', green)


%plot boundaries
% plot([0 4], [0 2], 'k')
% plot([0 -4], [0 2], 'k')

%% linearized version

%in non-parametric form, primary electrode axis = x, secondary electrode axis = f(x) = (x/lam) * ln(x/tp)
%first 2 terms of taylor series:
% f(x) = f(x_0) + f'(x_0)*(x-x_0) = (x_0/lam) * ln(x_0/tp) + (x-x_0)(1/lam)[ln(x_0/tp) + 1]
%
%

x_0_c1 = 0.35; %point to linearize around for cell 1 (in terms of primary electrode amplitude)
x_0_c2 = 0.58; %point to linearize around for cell 1
x_0_c3 = 0.78;

slope_1 = (1/lam_c1).*(log(x_0_c1/tp_c1)+1);
slope_2 = (1/lam_c2).*(log(x_0_c2/tp_c2)+1);
slope_3 = (1/lam_c3).*(log(x_0_c3/tp_c3)+1);

y0_1 = (x_0_c1/lam_c1)*log(x_0_c1/tp_c1) - x_0_c1*slope_1;
y0_2 = (x_0_c2/lam_c2)*log(x_0_c2/tp_c2) - x_0_c2*slope_2;
y0_3 = (x_0_c3/lam_c3)*log(x_0_c3/tp_c3) - x_0_c3*slope_3;

x = 0:0.01:2;

% threshold as a function of secondary:primary ratio for secondary electrode 2

t_1_fx = y0_1 + x.*slope_1;
t_2_fx = y0_2 + x.*slope_2;
t_3_fx = y0_3 + x.*slope_3;

plot(t_1_fx, x, 'LineWidth', 2, 'Color', blue)
plot(t_2_fx, x, 'LineWidth', 2, 'Color', green)
plot(x, t_3_fx, 'LineWidth', 2, 'Color', red) %axes switched because primary/secondary electrodes swapped

%plot standard deviations
% plot(t_1_fx*(1-1/w1), x*(1-1/w1), 'Color', blue)
% plot(t_1_fx*(1+1/w1), x*(1+1/w1), 'Color', blue)
% 
% plot(t_2_fx*(1-1/w2), x*(1-1/w2), 'Color', green)
% plot(t_2_fx*(1+1/w2), x*(1+1/w2), 'Color', green)
%  
% plot(x*(1-1/w3), t_3_fx*(1-1/w3), 'Color', red)
% plot(x*(1+1/w3), t_3_fx*(1+1/w3), 'Color', red)

%fill standard deviation regions

fillPolyY_1 = [xPlotLim(1) xPlotLim(2) xPlotLim(2) xPlotLim(1)];

fillPolyX_1 = zeros(4, 1);
fillPolyX_1(1:2) = (fillPolyY_1(1:2) - y0_1*(1+1/w1))/slope_1;
fillPolyX_1(3:4) = (fillPolyY_1(3:4) - y0_1*(1-1/w1))/slope_1;


fillPolyY_2 = [xPlotLim(1) xPlotLim(2) xPlotLim(2) xPlotLim(1)];

fillPolyX_2 = zeros(4, 1);
fillPolyX_2(1:2) = (fillPolyY_2(1:2) - y0_2*(1+1/w2))/slope_2;
fillPolyX_2(3:4) = (fillPolyY_2(3:4) - y0_2*(1-1/w2))/slope_2;


fillPolyY_3 = [xPlotLim(1) xPlotLim(2) xPlotLim(2) xPlotLim(1)];

fillPolyX_3 = zeros(4, 1);
fillPolyX_3(1:2) = (fillPolyY_3(1:2) - y0_3*(1+1/w3))/slope_3;
fillPolyX_3(3:4) = (fillPolyY_3(3:4) - y0_3*(1-1/w3))/slope_3;



fill(fillPolyY_1, fillPolyX_1, blue)
fill(fillPolyY_2, fillPolyX_2, green)
fill(fillPolyX_3, fillPolyY_3, red)


%locating point with equal stimulation on both electrodes on upper stdev bound for cell 1
bestInd = find(x*(1+1/w1) > 0.27, 1);
eqInd = find(x*(1+1/w1) > t_1_fx*(1+1/w1), 1);
plot(t_1_fx(eqInd)*(1+1/w1), x(eqInd)*(1+1/w1), 'ko', 'MarkerFaceColor', [0 0 0])
plot(t_1_fx(bestInd)*(1+1/w1), x(bestInd)*(1+1/w1), 'ko')



hold off
xlabel('electrode 1 current')
ylabel('electrode 2 current')
axis equal
set(gca, 'xlim', xPlotLim, 'ylim', yPlotLim, 'xtick', [0 0.5 1], 'ytick', [0 0.5 1])
