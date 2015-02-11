clear f_temp
syms p1 p2 tau1 tau2 t n

f = p1*((t/tau1)^n)*exp(-n*(t/tau1-1)) + p2*((t/tau2)^n)*exp(-n*(t/tau2-1))
%f = subs(f, p1, 0.1941);

f_diff = diff(f, t);

f_n = subs(f, n, 10);
f_diff_n = subs(f_diff, n, 10);


f_temp1 = subs(f_diff_n, t, 3);

f_temp2 = subs(f_diff_n, t, 7);

f_temp3 = subs(f_n, t, 4.75);

f_temp4 = subs(f_n, t, 3);
%f_temp3 = subs(f, t, 4.75);

[p1_init, p2_init, tau1_init, tau2_init] = solve(f_temp1==0, f_temp2 ==0, f_temp3 == 0, f_temp4 ==1, p1, p2, tau1, tau2)

f_n = subs(f_n, p1, p1_init);
f_n = subs(f_n, p2, p2_init);
f_n = subs(f_n, tau1, tau1_init);
f_n = subs(f_n, tau2, tau2_init);

f_plot = double(subs(f_n, t, [0:29]));
time = 30-[0:29];
figure
plot(time, f_plot)



