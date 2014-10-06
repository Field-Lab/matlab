%linearity toy model
clear all

% note: parameters have different form than erfFitter
b=2; %log(threshold)
w = 1;
x = 0:0.1:5;
R = 0.5*(1 + erf(w*(x-b)));

R2 = 0.5*(1 + erf(0.5*w*(x-2*b)));


%% constructing the toy model

lam1_c1 = -0.2; %cell one, secondary electrode 1
lam1_c2 = 0.2; %cell two, secondary electrode 1
lam2_c1 = -0.2; %cell one, secondary electrode 2
lam2_c2 = -0.1; %cell two, secondary electrode 2
tp_c1 = 0.5;
tp_c2 = 0.6;
w1 = 5;
w2 = 2.5;

c = -2:0.1:2; %secondary:primary ratios

blue = [27 117 187]/255;
red = [190 30 45]/255;

%% primary electrode + secondary electrode 1
figure %threshold as a function of secondary:primary ratio for secondary electrode 2
hold on
t_new_c1 = tp_c1*exp(c*lam1_c1); %cell 1
t_new_c2 = tp_c2*exp(c*lam1_c2); %cell 2
plot(t_new_c1, c.*t_new_c1, 'LineWidth', 2, 'Color', red)
plot(t_new_c2, c.*t_new_c2, 'LineWidth', 2, 'Color', blue)

%plot standard deviations
t_new_c1_sd = (t_new_c1/w1);
t_new_c2_sd = (t_new_c2/w2);

plot(t_new_c1 - t_new_c1_sd, c.*t_new_c1 - c.*t_new_c1_sd, 'Color', red)
plot(t_new_c1 + t_new_c1_sd, c.*t_new_c1 + c.*t_new_c1_sd, 'Color', red)

plot(t_new_c2 - t_new_c2_sd, c.*t_new_c2 - c.*t_new_c2_sd, 'Color', blue)
plot(t_new_c2 + t_new_c2_sd, c.*t_new_c2 + c.*t_new_c2_sd, 'Color', blue)

%plot boundaries
plot([0 2], [0 4], 'k')
plot([0 2], [0 -4], 'k')

hold off
xlabel('e_p')
ylabel('e_{s_1}')
axis equal
set(gca, 'xlim', [0 2], 'ylim', [-2 3])


%% primary electrode + secondary electrode 2
figure %threshold as a function of secondary:primary ratio for secondary electrode 2
hold on
t_new_c1 = tp_c1*exp(c*lam2_c1); %cell 1
t_new_c2 = tp_c2*exp(c*lam2_c2); %cell 2
plot(t_new_c1, c.*t_new_c1, 'LineWidth', 2, 'Color', red)
plot(t_new_c2, c.*t_new_c2, 'LineWidth', 2, 'Color', blue)

%plot standard deviations
t_new_c1_sd = (t_new_c1/w1);
t_new_c2_sd = (t_new_c2/w2);

plot(t_new_c1 - t_new_c1_sd, c.*t_new_c1 - c.*t_new_c1_sd, 'Color', red)
plot(t_new_c1 + t_new_c1_sd, c.*t_new_c1 + c.*t_new_c1_sd, 'Color', red)

plot(t_new_c2 - t_new_c2_sd, c.*t_new_c2 - c.*t_new_c2_sd, 'Color', blue)
plot(t_new_c2 + t_new_c2_sd, c.*t_new_c2 + c.*t_new_c2_sd, 'Color', blue)

%plot boundaries
plot([0 2], [0 4], 'k')
plot([0 2], [0 -4], 'k')

hold off
xlabel('e_p')
ylabel('e_{s_2}')
axis equal
set(gca, 'xlim', [0 2], 'ylim', [-2 2])





