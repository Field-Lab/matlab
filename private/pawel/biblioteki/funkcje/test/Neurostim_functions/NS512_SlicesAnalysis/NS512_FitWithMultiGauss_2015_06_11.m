function [CFs,leg,moje]=NS512_FitWithMultiGauss_2015_06_11(t,p);

fo = fitoptions('method','NonlinearLeastSquares','Normalize','off');
fo = fitoptions('method','NonlinearLeastSquares','Normalize','off','MaxFunEvals',10000,'MaxIter',10000);
ft = fittype('A*(exp(-(x-tau).^2/(2*sigma^2)))',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'A', 'tau', 'sigma'});
CFs=[];
window=ones(1,10); % dla 2015-05-24: window=ones(1,10)
for i=1:10            
    [c,smax,sig0]=ConvPeakMax(p,window);
    maxs=find(p==max(p));
    tau0=t(maxs(1));
    tau0=smax;
    sig0;
    %[max(p) tau0 2]
    set(fo,'Startpoint',[max(p) tau0 sig0]); %set(fo,'Startpoint',[max(p) tau0 2]);    
    [cf,g1] = fit(t',p',ft,fo);
    cf;
    if cf.sigma<0
        if i==1
            CFs=[];
            leg=[];
            moje=[];
        end
        return;
    end
    samples=round([floor(cf.tau-2*cf.sigma):ceil(cf.tau+2*cf.sigma)]);
    if (min(samples)<1 || max(samples)>600 || cf.sigma>100 || cf.A<0)
        'sdfsdfsdfsdf'
        if i==1
            CFs=[];
            leg=[];
            moje=[];
        end
        return;
    else
        'OK'
    end    
    
    cf(samples);
    %p(samples)
    p(samples)-cf(samples)';
    %(p(samples)-cf(samples)').^2
    %(p(samples)-cf(samples)').^2)./cf(samples)'
    chi2=sum(((p(samples)-cf(samples)').^2)./cf(samples)')/length(samples);
    %sum(cf(samples))
    moje(i)=sqrt(sum((p(samples)-cf(samples)').^2))/sum(cf(samples));
    Area=sqrt(6.28)*cf.A*(0.1+cf.sigma); % it is necessary to add 1 to sigma in case when the distribution is very narrow (for example, all spikes have identical delay)    
    if Area<50
        if i==1
            CFs=[];
            leg=[];
        end
        return
    end
    if i>1
        for j=1:i-1                  
            tau_old=CFs{j}.tau;
            sigma_old=CFs{j}.sigma;
            if (abs(tau_old-cf.tau))<3*(sigma_old+cf.sigma)
                cf.A=-1;            
            end
        end
    end
        
    CFs{i}=cf;    
    leg{i}=['t=' num2str(cf.tau,3) ', A=' num2str(Area,'%4.0f') ', m=' num2str(moje(i),3)];
    %leg{i}=['t=' num2str(cf.tau,3) ', A=' num2str(Area,'%4.0f') ', chi2=' num2str(chi2,3) ', moje=' num2str(moje(i),3)];
    t_min=max(round(cf.tau-3*cf.sigma),t(2)); %tutaj czasem daje zbyt duze wartosci, trzeba dodac jakis warunek typu "if"
    t_max=min(round(cf.tau+3*cf.sigma),t(length(t)));
    p([t_min:t_max])=0;    
    
    %TimeLocking=NS512_AreTheSpikesTimeLocked(p,window);
    %if TimeLocking==0
    %    return
    %end
    
    if find(p>0)
        
    else
        return;
    end
    %window=ones(1,40);
    %[c,smax]=ConvPeakMax(p,window)
    %pconv=conv(p,window);
    %conv_fin=pconv(length(window):length(p))/length(window);
    
    %if max(c)<3*mean(p)
    %    return;
    %end
end