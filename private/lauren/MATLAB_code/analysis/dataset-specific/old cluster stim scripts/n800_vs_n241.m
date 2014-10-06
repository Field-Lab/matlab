 cell_241_lam = [0.2949; %listed from top electrode, clockwise
                0.3543;
                0.2189;
                0.3434;
                0.1045;
                0.0950];
           
cell_800_lam = [0.4336;
               0.0878;
               0.1002;
               0.1407;
               -0.1757;
               0.1062];

           
          
          
pAloneThresh_241 = 2.1809/2;
pAloneThresh_800 = 0.4142;
          
w_241 = 2.1809;
w_800 = 3.4184;

secondaryAmpsSign = [-1 1 1 1 1 -1];


alpha = linspace(-1, 1, 10);

t_241_new = ones(size(alpha))*pAloneThresh_241;
t_800_new = ones(size(alpha))*pAloneThresh_800;
for i = 1:6
   t_241_new = t_241_new - cell_241_lam(i)*secondaryAmpsSign(i)*alpha;
   t_800_new = t_800_new - cell_800_lam(i)*secondaryAmpsSign(i)*alpha;
end

fillPolyX_241 = [alpha alpha(length(alpha):-1:1)];
fillPolyY_241 = [t_241_new*(1+1/w_241) t_241_new(length(alpha):-1:1)*(1-1/w_241)];

fillPolyX_800 = [alpha alpha(length(alpha):-1:1)];
fillPolyY_800 = [t_800_new*(1+1/w_241) t_800_new(length(alpha):-1:1)*(1-1/w_800)];

figure
hold on
fill(fillPolyX_241, fillPolyY_241, [1 0.8 0.8])
fill(fillPolyX_800, fillPolyY_800, [0.8 0.8 1])
plot(alpha, t_241_new, 'LineWidth', 1, 'Color', 'r')
plot(alpha, t_800_new, 'LineWidth', 1, 'Color', 'b')
hold off
axis equal
set(gca, 'xlim', [-1 1], 'ylim', [0 2])
