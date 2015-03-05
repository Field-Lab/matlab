%this function defines conditions of differential equations, calls
%computational function and creates voltage plots
global Vth
Vth=-50; %mV
tspan = [0 200];
y0 = -65;
% Solve the problem using ode45 (based on Runge-Kutta, medium accuracy)
solution=ode45(@myIAF,tspan,y0);
%min and max values for y axis
ymin=-70;
ymax=-20;
%plot solution of chosen compartment as potential vs time, and add stimulus mark to appropriate
%plot.
figure(2)
plot(solution.x,solution.y')
axis([0 200 ymin ymax])
