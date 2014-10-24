clf
a = 0.1;
b = 1/a;
c = 0;
d = 1;
e = 0;

% when c is 0, d is 1, a is small, and b is 1/a
% this is approximately rectified threshold

% when c is 0, d is 0, and a is 1, it's linear
%%
clf
x = [-5:0.01:5];
c = 4;
b = 0.1;
a = 0;
%y = c*normcdf(b*x + a);
y = c * log(1+exp(b*x+a));
plot(x,y);
xlim([-5 5]);
ylim([-5 5]);
drawHorzLine(0);
drawVertLine(0);

%%


%%
clf
x = [-5:0.01:10];
f.type = 'loglinear';
f.p = [0.4 3 0.1 0];
    
[y df dff] = evalNonLin(x,f);

plot(x,y);
hold on
plot(x,df,'r');
plot(x,dff,'k');
xlim([-5 10]);
ylim([-5 40]);
drawHorzLine(0);
drawVertLine(0);

hold on
%plot(x,max(x,0).^2,'r');
%plot(x,max(x,0),'g');

%%
clf
x = [-5:0.01:5];
sim.type = 'exp';
sim.p = [0 1];
y = evalNonLin(x,sim);
%y = y + randn(size(y))*0.2;
nl.type = 'polynomial';
nl = fitNLpar(x,y,nl);
nl.p
%nl.p = [0 1 10 0.1 -12];
hold on
plot(x,y,'k','LineWidth',10);
plot(x,evalNonLin(x,nl),'r','LineWidth',3);
ylim([-5 5]);
xlim([-5 5]);
%%
x = [-5:0.01:5];
y1 = x;
y2 = x.^2;
y3 = x.^3;
a = 15;
b = 1;
c = 0.5;
y = a*y1 + b*y2 + c*y3;

plot(x,y);

%%
cellTypes = {'off midget'};
iType = 1;
clf
set(gcf,'Position',[88 114 1299 692]);
ids = getTypeIds(dat.cellTypes,cellTypes{iType});
r2_SUB = zeros(1,length(ids));
r2_LN = zeros(1,length(ids));

%%
clf
%ids = shuffle(ids);
plotIter = 1;
i = 1;

while plotIter <= 40 && i <= length(ids)
    
    try
        fileName = sprintf('%g-subunit-cliques-20-fit.mat',ids(i));
        if isfile(fullfile(dat.analPath,'subunit',fileName));
            load(fullfile(dat.analPath,'subunit',fileName));
        else
            fileName = sprintf('%g-subunit-singletons-fit.mat',ids(i));
            load(fullfile(dat.analPath,'subunit',fileName));
        end
        r2_SUB(i) = res.out_SUB.r2;
        r2_LN(i) = res.out_LN.r2;
        if r2_SUB(i) > 0.05;
            subp(5,8,plotIter,0.02);
            nl_spline = res.fit_SUB.f;
            x = linspace(min(nl_spline.knots),max(nl_spline.knots),100);
            y = evalNonLin(x,nl_spline);
            nl.type = 'polynomial';
            nl = fitNLpar(x,y,nl);
            hold on
            yHat = evalNonLin(x,nl);
            sse = sum((yHat-y).^2);
            sst = sum((y-mean(y)).^2);
            r2(plotIter) = 1 - sse / sst
            plot(x,yHat,'Color',[0.2 0.7 0.7],'LineWidth',5);
            plotNonLin(nl_spline);
            plotIter = plotIter + 1;
            xlim([-2.5 2.5]);
            ylim([-1 2.75]);
            drawnow;
        end
        
    catch
        fprintf('(batchNonLin) cannot load cell %g\r',i);
    end
    i = i + 1;
    drawnow;
end