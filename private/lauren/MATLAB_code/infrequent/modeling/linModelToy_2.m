%linearity toy model 2: 2 cells between 2 electrodes, with 2 additional flanking electrodes

clear all

%% constructing the toy model

lam1_c1 = -0.2; %cell one, secondary electrode
lam1_c2 = -0.3; %cell two, secondary electrode

lamp_c1 = -0.5; %cell one, when secondary electrode acts as primary electrode
lamp_c2 = -0.4;

tp_c1 = 0.5;
ts_c1 = 0.9; %cell 1, threshold to secondary electrode
tp_c2 = 0.6;
ts_c2 = 1.2;

w1 = 5;
w2 = 3;

w1s = 4;

w2s = 2;

blue = [27 117 187]/255;
red = [190 30 45]/255;

%% primary electrode + secondary electrode 1
figure %threshold as a function of secondary:primary ratio for secondary electrode 2
hold on

%cell 1: region where primary acts as primary
c = -2:0.1:ts_c1/tp_c1;
t_new_c1 = tp_c1*exp(c*lam1_c1); %cell 1
plot(t_new_c1, c.*t_new_c1, 'LineWidth', 2, 'Color', red)

%standard deviations
t_new_c1_sd = (t_new_c1/w1);
plot(t_new_c1 - t_new_c1_sd, c.*t_new_c1 - c.*t_new_c1_sd, 'Color', red)
plot(t_new_c1 + t_new_c1_sd, c.*t_new_c1 + c.*t_new_c1_sd, 'Color', red)

%cell 1: region where secondary acts as primary (lower boundary for c is unknown)
c = -1:0.1:tp_c1/ts_c1;
t_new_c1 = ts_c1*exp(c*lamp_c1); %cell 1
plot(c.*t_new_c1, t_new_c1, 'LineWidth', 2, 'Color', 0.7*red)

%standard deviations
t_new_c1_sd = (t_new_c1/w1s);
plot(c.*t_new_c1 - c.*t_new_c1_sd, t_new_c1 - t_new_c1_sd, 'Color', 0.7*red)
plot(c.*t_new_c1 + c.*t_new_c1_sd, t_new_c1 + t_new_c1_sd, 'Color', 0.7*red)

%plot boundaries
plot([0 2], [0 2*ts_c1/tp_c1], 'Color', red)
plot([0 2], [0 -4], 'Color', red)
plot([0 -2], [0 2], '--','Color', 0.7*red)


%cell 2: axes are flipped because primary/secondary are switched
c = -2:0.1:ts_c2/tp_c2;
t_new_c2 = tp_c2*exp(c*lam1_c2); %cell 2
plot(c.*t_new_c2, t_new_c2, 'LineWidth', 2, 'Color', blue)

%standard deviations
t_new_c2_sd = (t_new_c2/w2);
plot(c.*t_new_c2 - c.*t_new_c2_sd, t_new_c2 - t_new_c2_sd, 'Color', blue)
plot(c.*t_new_c2 + c.*t_new_c2_sd, t_new_c2 + t_new_c2_sd, 'Color', blue)

%cell 2: region where secondary acts as primary
c = -1:0.1:tp_c2/ts_c2;
t_new_c2 = ts_c2*exp(c*lamp_c2); %cell 2
plot(t_new_c2, c.*t_new_c2, 'LineWidth', 2, 'Color', 0.7*blue)

%standard deviations
t_new_c2_sd = (t_new_c2/w2s);
plot(t_new_c2 - t_new_c2_sd, c.*t_new_c2 - c.*t_new_c2_sd, 'Color', 0.7*blue)
plot(t_new_c2 + t_new_c2_sd, c.*t_new_c2 + c.*t_new_c2_sd, 'Color', 0.7*blue)

%plot boundariescd 
plot([0 2*ts_c2/tp_c2], [0 2], 'Color', blue)
plot([0 -4], [0 2], 'Color', blue)
plot([0 2], [0 -2], '--','Color', 0.7*blue)

hold off
xlabel('e_p (cell 1), e_s (cell 2)')
ylabel('e_s (cell 1), e_p (cell 2)')
axis equal
set(gca, 'xlim', [-3 3], 'ylim', [-2 3])