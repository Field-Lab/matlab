function [pred,lam,nl_info]=response_lnp_optimal_nl(data,stimulus,nTrials,recorded_resp)

 mov=stimulus; % 3D, between -0.5 to +0.5.
maskedMovdd= filterMov(mov,data.mask,squeeze(data.ttf));
inputt = (data.usta'*maskedMovdd - data.ainput) / data.binput;


% phase 1
input = inputt'; output = recorded_resp'/(1/120);
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
[fitres1, gof] = fit(inputt',recorded_resp'/(1/120), 'scale*normcdf(x,mu,sigma)', 'StartPoint', [fitres.mu, fitres.scale, fitres.sigma]);

figure;
plot(inputt',recorded_resp'/(1/120),'.');
hold on;
plot(in',fitres.scale*normcdf(in,fitres.mu,fitres.sigma),'LineWidth',3);
hold on;
plot(in',fitres1.scale*normcdf(in,fitres1.mu,fitres1.sigma),'LineWidth',3);

lam = fitres1.scale*normcdf(inputt,fitres1.mu,fitres1.sigma);
lam = lam*data.dt;


pred=zeros(nTrials,length(lam));
for itrial=1:nTrials
pred(itrial,:)=poissrnd(lam);    
end

nl_info.input = inputt';
nl_info.recorded_resp = recorded_resp'/(1/120);
nl_info.fit_params = fitres1;
nl_info.nl_expression = 'scale*normcdf(x,mu,sigma)'; 

end