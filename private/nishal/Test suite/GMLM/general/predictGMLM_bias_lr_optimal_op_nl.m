
function [pred,lam,nl_info] = predictGMLM_bias_lr_optimal_op_nl(fitGMLM,maskedMov,nTrials,interval,recorded_resp)
dt = 1/120;
[~,lam]=predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,interval);

% phase 1
input = lam*120; output = recorded_resp'/(1/120);
thrr=[];
for iprc=0:2:100
thrr=[thrr;prctile(input,iprc)];
end

fr = [];
in =[];
err=[];
for ipt = 1:length(thrr)-1
    idx = input>thrr(ipt) & input<=thrr(ipt+1);
    fr = [fr;mean(output(idx))];
    err = [err;sqrt(var(output(idx))/sum(idx))];
    in = [in;mean(input(idx))];
end

fitres = fit(in,fr, 'scale*normcdf(x,mu,sigma)','StartPoint',[1, max(fr), 1.5]);
[fitres1, gof] = fit(input,output, 'scale*normcdf(x,mu,sigma)', 'StartPoint', [fitres.mu, fitres.scale, fitres.sigma]);

figure;
semilogx(input,output,'.');
hold on;
semilogx(in',fitres.scale*normcdf(in,fitres.mu,fitres.sigma),'LineWidth',3);
hold on;
semilogx(in',fitres1.scale*normcdf(in,fitres1.mu,fitres1.sigma),'LineWidth',3);

lam = fitres1.scale*normcdf(input,fitres1.mu,fitres1.sigma);
lam = lam/120;


pred=zeros(nTrials,length(lam));
for itrial=1:nTrials
pred(itrial,:)=poissrnd(lam);    
end

nl_info.input = input;
nl_info.recorded_resp =output;
nl_info.fit_params = fitres1;
nl_info.nl_expression = 'scale*normcdf(x,mu,sigma)'; 

end