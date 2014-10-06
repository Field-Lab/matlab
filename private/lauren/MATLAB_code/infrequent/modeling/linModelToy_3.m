%linearity toy model 2: 2 cells between 2 electrodes, with 2 additional flanking electrodes

clear all

%% constructing the toy model

lam1_c1 = 0.05; %cell one, secondary electrode
lam1_c2 = -0.1; %cell two, secondary electrode
lam1_c3 = -0.25;

tp_c1 = 0.5;
tp_c2 = 0.7;
tp_c3 = 0.6;

w1 = 5;
w2 = 3;
w3 = 2;

blue = [27 117 187]/255;
red = [190 30 45]/255;
green = [41 180 115]/255;

%% primary electrode + secondary electrode 1
figure %threshold as a function of secondary:primary ratio for secondary electrode 2
hold on

%cell 1
c = -2:0.1:2;
t_new_c1 = tp_c1*exp(c*lam1_c1); %cell 1
plot(t_new_c1, c.*t_new_c1, 'LineWidth', 2, 'Color', red)

%standard deviations
t_new_c1_sd = (t_new_c1/w1);
plot(t_new_c1 - t_new_c1_sd, c.*t_new_c1 - c.*t_new_c1_sd, 'Color', red)
plot(t_new_c1 + t_new_c1_sd, c.*t_new_c1 + c.*t_new_c1_sd, 'Color', red)


%plot boundaries
plot([0 2], [0 4], 'Color', red)
plot([0 2], [0 -4], 'Color', red)



%cell 2: axes are flipped because primary/secondary are switched
c = -2:0.1:2;
t_new_c2 = tp_c2*exp(c*lam1_c2); %cell 2
plot(c.*t_new_c2, t_new_c2, 'LineWidth', 2, 'Color', blue)

%standard deviations
t_new_c2_sd = (t_new_c2/w2);
plot(c.*t_new_c2 - c.*t_new_c2_sd, t_new_c2 - t_new_c2_sd, 'Color', blue)
plot(c.*t_new_c2 + c.*t_new_c2_sd, t_new_c2 + t_new_c2_sd, 'Color', blue)

%plot boundaries
plot([0 4], [0 2], 'Color', blue)
plot([0 -4], [0 2], 'Color', blue)



%cell 3
c = -2:0.1:2;
t_new_c3 = tp_c3*exp(c*lam1_c3); %cell 1
plot(t_new_c3, c.*t_new_c3, 'LineWidth', 2, 'Color', green)

%standard deviations
t_new_c3_sd = (t_new_c3/w3);
plot(t_new_c3 - t_new_c3_sd, c.*t_new_c3 - c.*t_new_c3_sd, 'Color', green)
plot(t_new_c3 + t_new_c3_sd, c.*t_new_c3 + c.*t_new_c3_sd, 'Color', green)

%plot boundaries
plot([0 2], [0 4], '--', 'Color', green)
plot([0 2], [0 -4], '--', 'Color', green)



hold off
xlabel('e_p (cells 1,3), e_s (cell 2)')
ylabel('e_s (cells 1,3), e_p (cell 2)')
axis equal
set(gca, 'xlim', [-2 2], 'ylim', [-2 2])

