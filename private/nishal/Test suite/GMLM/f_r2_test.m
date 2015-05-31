function [f_val_test,R2_log]=f_r2_test(lam,rec,interval)
       %% R2 value 
%       predd=lam; 
%      
%        pred_ss = zeros(length(predd)/10,1);
%        for itime=1:length(pred_ss)
%        pred_ss(itime) = sum(predd((itime-1)*10+1:(itime)*10));
%        end
       pred_ss = lam';
       
       % R2 value method 2
       x1 = pred_ss; y1 = rec; n=length(y1);
       r = (n*x1'*y1 - sum(x1)*sum(y1))/(sqrt(n*sum(x1.^2) - sum(x1)^2) * sqrt(n*sum(y1.^2) - sum(y1)^2));
       R2_log = r^2;
       
       
       %% Likelihood value
       lam = pred_ss;
       tsp = find(rec~=0);
       y_tsp = rec(rec~=0);

        lam=lam;

       likelihood = (sum(-lam*(interval/1200)) + sum(y_tsp.*log(lam(tsp)*(interval/1200))))/length(rec);
       f_val_test=-likelihood;
end