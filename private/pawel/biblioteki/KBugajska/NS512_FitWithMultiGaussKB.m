function CFs=NS512_FitWithMultiGauss(t,p);

fo = fitoptions('method','NonlinearLeastSquares','Normalize','off');
fo = fitoptions('method','NonlinearLeastSquares','Normalize','off','MaxFunEvals',10000,'MaxIter',10000);
ft = fittype('A*(exp(-(x-tau).^2/(2*sigma^2)))',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'A', 'tau', 'sigma'});
CFs=[];
for i=1:10
    maxs=find(p==max(p));
    tau0=t(maxs(1));
    [max(p) tau0 2]
    set(fo,'Startpoint',[max(p) tau0 2]); %set(fo,'Startpoint',[max(p) tau0 2]);
    cf = fit(t',p',ft,fo);
    Area=sqrt(6.28)*cf.A*(1+cf.sigma); % it is necessary to add 1 to sigma in case when the distribution is very narrow (for example, all spikes have identical delay)    
    if i>1
        for j=1:i-1       
            if cf.A<0
                return;
            end
            if Area<100
                return;
            end
            tau_old=CFs{j}.tau;
            sigma_old=CFs{j}.sigma;
            if (abs(tau_old-cf.tau))<1.1*(sigma_old+cf.sigma)
                return;
            end
        end
    end    
    
    CFs{i}=cf;
    
   % p=NS512_VanishingDataBySigmaKB(t,p,cf);
  
   switch logical(true)
            
        case (cf.sigma<1)
            
          t_min=max(round(cf.tau-(20*cf.sigma)),1);
          t_max=min(round(cf.tau+(20*cf.sigma)), t(length(t)));
        
           p(t_min:t_max)=0;
             
        case cf.sigma<2  
            
          t_min=max(round(cf.tau-(12*cf.sigma)),1);
          t_max=min(round(cf.tau+(12*cf.sigma)),  t(length(t)));
        
           p(t_min:t_max)=0;
         
        case cf.sigma<10
            
         t_min=max(round(cf.tau-(10*cf.sigma)),1);
         t_max=min(round(cf.tau+(10*cf.sigma)),  t(length(t)));
        
            p(t_min:t_max)=0;
        
        case cf.sigma <50
            
         t_min=max(round(cf.tau-(6*cf.sigma)),1);
         t_max=min(round(cf.tau+(6*cf.sigma)), t(length(t)));
        
            p(t_min:t_max)=0;
        
        case cf.sigma <100 
            
         t_min=max(round(cf.tau-(4*cf.sigma)),1);
         t_max=min(round(cf.tau+(4*cf.sigma)), t(length(t)));
        
          p(t_min:t_max)=0;
        
        case cf.sigma >=100 
            
         t_min=max(round(cf.tau-(3*cf.sigma)),1);
         t_max=min(round(cf.tau+(3*cf.sigma)), t(length(t)));
        
         p(t_min:t_max)=0;
      
        otherwise warning('Unexpected sigma');

    end
    
end
   


