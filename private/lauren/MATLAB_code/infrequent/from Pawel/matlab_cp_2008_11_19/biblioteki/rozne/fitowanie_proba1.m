x = (0:100)';
y = 5 + 10./(1+exp(-(x-40)/10)) + randn(size(x));
plot(x,y,'bo')

%%
f = @(p,x) p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4)));
p = nlinfit(x,y,f,[10 20 50 5])
line(x,f(p,x),'color','r')

%%
fit(x,y,'a + b ./ (1 + exp(-(x-m)/s))','start',[0 20 50 5]) 