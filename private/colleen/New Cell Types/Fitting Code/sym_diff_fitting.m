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
[W,X,Y,Z] = ndgrid([0.1:.05:0.3],[-0.1:.01:-0.01],[1:1:10], [1:1:10]);

f_temp1_eval = zeros(size(W));
f_temp2_eval = zeros(size(W));
f_temp3_eval = zeros(size(W));
f_temp4_eval = zeros(size(W));
for h = 1:size(X,1)
    for i = 1:size(X,2)
        for j = 1:size(X,3)
            for k = 1:size(X,4)
                f_temp1_eval(h,i,j,k) = subs(f_temp1,{p1,p2, tau1, tau2},[W(h,i,j,k), X(h,i,j,k),Y(h,i,j,k),Z(h,i,j,k)]);
                f_temp2_eval(h,i,j,k) = subs(f_temp2,{p1,p2, tau1, tau2},[W(h,i,j,k), X(h,i,j,k),Y(h,i,j,k),Z(h,i,j,k)]);
                f_temp3_eval(h,i,j,k) = subs(f_temp3,{p1,p2, tau1, tau2},[W(h,i,j,k), X(h,i,j,k),Y(h,i,j,k),Z(h,i,j,k)]);
                f_temp4_eval(h,i,j,k) = subs(f_temp4,{p1,p2, tau1, tau2},[W(h,i,j,k), X(h,i,j,k),Y(h,i,j,k),Z(h,i,j,k)])-1;

            end
        end
    end
    h
end
temp1_norm = (f_temp1_eval - min(f_temp1_eval(:))) / ( max(f_temp1_eval(:)) - min(f_temp1_eval(:)) );
temp2_norm = (f_temp2_eval - min(f_temp2_eval(:))) / ( max(f_temp2_eval(:)) - min(f_temp2_eval(:)) );
temp3_norm = (f_temp3_eval - min(f_temp3_eval(:))) / ( max(f_temp3_eval(:)) - min(f_temp3_eval(:)) );
temp4_norm = (f_temp4_eval - min(f_temp4_eval(:))) / ( max(f_temp4_eval(:)) - min(f_temp4_eval(:)) );

%f_temp4_eval(f_temp4_eval <0) = 100;
temp1_norm(temp1_norm ==0 ) = 100;
temp2_norm(temp2_norm ==0 ) = 100;
temp3_norm(temp3_norm ==0 ) = 100;

temp_norm = temp1_norm + temp2_norm + temp3_norm + temp4_norm;
%f_temp_eval = f_temp1_eval + f_temp2_eval + f_temp3_eval + f_temp4_eval;
[x,y] = min(abs(temp_norm(:)))

[a,b,c,d] = ind2sub(size(W), y)
p1_init = W(a,b,c,d);
p2_init = X(a,b,c,d);
tau1_init = Y(a,b,c,d);
tau2_init = Z(a,b,c,d);
f_temp2 = subs(f_temp2, p1, 1);


f_temp3 = subs(f_n, t, 4.75);
f_temp3 = subs(f_temp3, p1, 1);

%f_temp4 = subs(f_n, t, 3);





%f_temp3 = subs(f, t, 4.75);

[tau1_init] = solve(f_temp1<=0, tau1)

f_n = subs(f_n, p1, p1_init);
f_n = subs(f_n, p2, p2_init);
f_n = subs(f_n, tau1, tau1_init);
f_n = subs(f_n, tau2, tau2_init);

f_plot = double(subs(f_n, t, [0:29]));
time = 30-[0:29];
figure
plot(time, f_plot)



