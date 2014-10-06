%linearity toy model
clear all

% note: parameters have different form than erfFitter
b=2; %log(threshold)
w = 1;
x = 0:0.1:5;
R = 0.5*(1 + erf(w*(x-b)));

R2 = 0.5*(1 + erf(0.5*w*(x-2*b)));


%% constructing the toy model

lam1_c1 = -0.3; %cell one, secondary electrode 1
lam1_c2 = 0; %cell two, secondary electrode 1
lam2_c1 = -0.05; %cell one, secondary electrode 2
lam2_c2 = -0.2; %cell two, secondary electrode 2
tp_c1 = 0.5;
tp_c2 = 0.6;
w1 = 5;
w2 = 2.5;

alpha = -2:0.1:2; %secondary:primary ratios

blue = [27 117 187]/255;
red = [190 30 45]/255;

%% primary electrode + secondary electrode 1
figure %threshold as a function of secondary:primary ratio for secondary electrode 2
hold on
t_new_c1 = tp_c1*exp(alpha*lam1_c1); %cell 1
t_new_c2 = tp_c2*exp(alpha*lam1_c2); %cell 2
plot(alpha.*t_new_c1, t_new_c1, 'LineWidth', 2, 'Color', red)
plot(alpha.*t_new_c2, t_new_c2, 'LineWidth', 2, 'Color', blue)

%plot standard deviations
t_new_c1_sd = (t_new_c1/w1);
t_new_c2_sd = (t_new_c2/w2);

plot(alpha.*t_new_c1*(1-1/w1), t_new_c1*(1-1/w1), 'Color', red)
plot(alpha.*t_new_c1*(1+1/w1), t_new_c1*(1+1/w1), 'Color', red)
%plot(t_new_c1 - t_new_c1_sd, alpha.*t_new_c1 - alpha.*t_new_c1_sd, 'Color', red)
%plot(t_new_c1 + t_new_c1_sd, alpha.*t_new_c1 + alpha.*t_new_c1_sd, 'Color', red)

plot(alpha.*t_new_c2*(1-1/w2), t_new_c2*(1-1/w2), 'Color', blue)
plot(alpha.*t_new_c2*(1+1/w2), t_new_c2*(1+1/w2), 'Color', blue)
%plot(t_new_c2 - t_new_c2_sd, alpha.*t_new_c2 - alpha.*t_new_c2_sd, 'Color', blue)
%plot(t_new_c2 + t_new_c2_sd, alpha.*t_new_c2 + alpha.*t_new_c2_sd, 'Color', blue)

%plot boundaries
plot([0 4], [0 2], 'k')
plot([0 -4], [0 2], 'k')

hold off
ylabel('i_p')
xlabel('i_{s_1}')
axis equal
set(gca, 'xlim', [-2 2], 'ylim', [0 2])


%% primary electrode + secondary electrode 2
figure %threshold as a function of secondary:primary ratio for secondary electrode 2
hold on
t_new_c1 = tp_c1*exp(alpha*lam2_c1); %cell 1
t_new_c2 = tp_c2*exp(alpha*lam2_c2); %cell 2
plot(alpha.*t_new_c1, t_new_c1, 'LineWidth', 2, 'Color', red)
plot(alpha.*t_new_c2, t_new_c2, 'LineWidth', 2, 'Color', blue)

%plot standard deviations
plot(alpha.*t_new_c1*(1-1/w1), t_new_c1*(1-1/w1), 'Color', red)
plot(alpha.*t_new_c1*(1+1/w1), t_new_c1*(1+1/w1), 'Color', red)
% plot(t_new_c1 - t_new_c1_sd, alpha.*t_new_c1 - alpha.*t_new_c1_sd, 'Color', red)
% plot(t_new_c1 + t_new_c1_sd, alpha.*t_new_c1 + alpha.*t_new_c1_sd, 'Color', red)

plot(alpha.*t_new_c2*(1-1/w2), t_new_c2*(1-1/w2), 'Color', blue)
plot(alpha.*t_new_c2*(1+1/w2), t_new_c2*(1+1/w2), 'Color', blue)
% plot(t_new_c2 - t_new_c2_sd, alpha.*t_new_c2 - alpha.*t_new_c2_sd, 'Color', blue)
% plot(t_new_c2 + t_new_c2_sd, alpha.*t_new_c2 + alpha.*t_new_c2_sd, 'Color', blue)

%plot boundaries
plot([0 4], [0 2], 'k')
plot([0 -4], [0 2], 'k')

hold off
ylabel('i_p')
xlabel('i_{s_2}')
axis equal
set(gca, 'xlim', [-2 2], 'ylim', [0 2])

%% primary electrode + both secondary electrodes
alpha1 = -2:0.2:2; %ratio of secondary 1 to primary
alpha2 = -2:0.2:2; %ratio of secondary 2 to primary

t_new_c1 = zeros(length(alpha1), length(alpha2));
t_new_c2 = zeros(length(alpha1), length(alpha2));
x1_c1 = zeros(length(alpha1), length(alpha2));
x1_c2 = zeros(length(alpha1), length(alpha2));
x2_c1 = zeros(length(alpha1), length(alpha2));
x2_c2 = zeros(length(alpha1), length(alpha2));

cValuesRed = zeros(length(alpha1), length(alpha2), 3);
cValuesBlue = zeros(length(alpha1), length(alpha2), 3);
for i = 1:length(alpha1)
    for j = 1:length(alpha2)
        cValuesRed(i,j,:) = red;
        cValuesBlue(i,j,:) = blue;        
    end
end

for i = 1:length(alpha2)
    %z-values (primary electrode dimension)
    t_new_c1(:,i) = tp_c1*exp(alpha1*lam1_c1 + alpha2(i)*lam2_c1); %cell 1
    t_new_c2(:,i) = tp_c2*exp(alpha1*lam1_c2 + alpha2(i)*lam2_c2); %cell 2
    
    %x1-values (secondary electrode 1 dimension)
    x1_c1(:,i) = alpha1'.*t_new_c1(:,i);
    x1_c2(:,i) = alpha1'.*t_new_c2(:,i);
    
    %x2-values (secondary electrode 2 dimenstion)
    x2_c1(:,i) = alpha2(i)*t_new_c1(:,i);
    x2_c2(:,i) = alpha2(i)*t_new_c2(:,i);
    
end

figure
hold on

mesh(x1_c1, x2_c1, t_new_c1, cValuesRed, 'FaceColor', 'none')
set(gco,'LineSmoothing','on') 
mesh(x1_c2, x2_c2, t_new_c2, cValuesBlue, 'FaceColor', 'none')
set(gco,'LineSmoothing','on') 

%standard deviations (scale by 1 - 1/w and 1 + 1/w)
mesh(x1_c1*(1-1/w1), x2_c1*(1-1/w1), t_new_c1*(1-1/w1), 1.5*cValuesRed, 'FaceColor', 'none')
set(gco,'LineSmoothing','on') 
mesh(x1_c1*(1+1/w1), x2_c1*(1+1/w1), t_new_c1*(1+1/w1), 0.5*cValuesRed, 'FaceColor', 'none')

mesh(x1_c2*(1-1/w2), x2_c2*(1-1/w2), t_new_c2*(1-1/w2), 1.5*cValuesBlue, 'FaceColor', 'none')
mesh(x1_c2*(1+1/w2), x2_c2*(1+1/w2), t_new_c2*(1+1/w2), 0.5*cValuesBlue, 'FaceColor', 'none')


plot3([0 0], [0 0], [0 1.2], 'k-')
plot3([0 0], [-2 2], [0 0], 'k-')
plot3([-2 2], [0 0], [0 0], 'k-')
plot3([0 3], [0 3], [0 1.5], 'r-')
plot3([0 -3], [0 3], [0 1.5], 'r-')
plot3([0 3], [0 -3], [0 1.5], 'r-')
plot3([0 -3], [0 -3], [0 1.5], 'r-')



set(gcf, 'color', 'white')




