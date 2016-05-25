function CFs=NS512_FitWithMultiGauss_2015_05_24(t,p);

fo = fitoptions('method','NonlinearLeastSquares','Normalize','off');
fo = fitoptions('method','NonlinearLeastSquares','Normalize','off','MaxFunEvals',10000,'MaxIter',10000);
ft = fittype('A*(exp(-(x-tau).^2/(2*sigma^2)))',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'A', 'tau', 'sigma'});
CFs=[];
for i=1:10        
    maxs=find(p==max(p));
    tau0=t(maxs(1));
    %[max(p) tau0 2]
    set(fo,'Startpoint',[max(p) tau0 2]); %set(fo,'Startpoint',[max(p) tau0 2]);
    %size(t)
    %size(p)
    cf = fit(t',p',ft,fo);
    Area=sqrt(6.28)*cf.A*(1+cf.sigma); % it is necessary to add 1 to sigma in case when the distribution is very narrow (for example, all spikes have identical delay)    
    if i>1
        for j=1:i-1       
            if cf.A<0
                return;
            end
            if Area<50 % was 100! 2013-10-10
                return;
            end
            %tau_new=cf.tau
            %sigma_new=cf.sigma
            tau_old=CFs{j}.tau;
            sigma_old=CFs{j}.sigma;
            if (abs(tau_old-cf.tau))<2*(sigma_old+cf.sigma)
                return;
            end
        end
    end
    
    CFs{i}=cf;
    t_min=max(round(cf.tau-3*cf.sigma),t(2)); %tutaj czasem daje zbyt duze wartosci, trzeba dodac jakis warunek typu "if"
    t_max=min(round(cf.tau+3*cf.sigma),t(length(t)));
    p([t_min:t_max])=0;    
    
    window=ones(1,40);
    [c,smax]=ConvPeakMax(p,window)
    %pconv=conv(p,window);
    %conv_fin=pconv(length(window):length(p))/length(window);
    
    if max(c)<3*mean(p)
        return;
    end
end