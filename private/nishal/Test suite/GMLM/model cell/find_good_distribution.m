
dy=0.01;

y_list=[dy:dy:1-dy];
mu =0;
sigma=2;



g=@(x) mu+sqrt(2)*sigma*erfinv(2*x- 1);
x_normal=g(y_list);
icnt=1;
for S=3:1:15
C = sum(exp(g(y_list))*dy);


iy1=0;
for  y1=dy:dy:1-dy
    iy1=iy1+1;
    iy2=0;
    for y2=dy:dy:1-dy
        iy2=iy2+1;
        p(iy1,iy2)=(exp(g(y1)) + exp(g(y2)) + (S-2)*C- 1.1*(S)) *((y2)^(S-2)) *(S*(S-1));
        valid(iy1,iy2)=(y1>y2);
    end
end

p=p.*valid;
p=p/sum(p(:));

y1_marginal =sum(p,2);
y1_marginal = y1_marginal / sum(y1_marginal);

y2_marginal =sum(p,1)';
y2_marginal = y2_marginal / sum(y2_marginal);

figure;
plot(y_list,y1_marginal,'r');
hold on;
plot(y_list,y2_marginal,'g');
legend('Max','Second max');

figure;
contour(y_list,y_list,p,20)


figure;
contour3(x_normal,x_normal,p,20);

meanact1 = exp(g(y_list))*y1_marginal;
meanact2 = exp(g(y_list))*y2_marginal;

diffmean(icnt)=meanact1-meanact2;
sCnt(icnt)=S;
icnt=icnt+1;

end

figure;
plot(sCnt,diffmean);