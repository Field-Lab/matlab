function sig_log = plot_cone_interaction(binnedResponsesbigd  , maskedMov ,cones,plot_fits)

icone = cones(1);jcone=cones(2);sig_log=[];
col='rrrrrrrrrrrrr'; col(jcone)='b';
%figure;
% subplot(3,1,1);

nbins=7;
for ibin=1:nbins+1
X(ibin) = prctile(maskedMov(icone,:),100*(ibin-1)/nbins);
end

spkRate=[];
for iconeInp=1:nbins-1
    idx = (maskedMov(icone,:) >= X(iconeInp) & maskedMov(icone,:) <= X(iconeInp+1));
spkRate(iconeInp) = sum(binnedResponsesbigd(idx))/sum(idx);
iconeInp_log(iconeInp) = mean(maskedMov(icone,idx));
end
fit_res = fit(iconeInp_log',spkRate','normcdf(x,mu,sigma)+b');

sig_log=[sig_log;fit_res.sigma];

if(~plot_fits)
plot(iconeInp_log,spkRate,'k');
hold on;
else
plot(iconeInp_log',feval(fit_res,iconeInp_log'),'k') 
hold on;
end

for jcone=[1:icone-1 ,icone+1:size(maskedMov,1)]
%     pause
thr = prctile(maskedMov(jcone,:),90);

for ibin=1:nbins+1
X(ibin) = prctile(maskedMov(icone, maskedMov(jcone,:)>=thr),100*(ibin-1)/nbins);
end


iconeInp_log=[];
spkRate=[];
for iconeInp=1:nbins-1
    
    idx = (maskedMov(icone,:) > X(iconeInp) & maskedMov(icone,:) <= X(iconeInp+1) & maskedMov(jcone,:)>=thr);
spkRate(iconeInp) = sum(binnedResponsesbigd(idx))/sum(idx);
iconeInp_log(iconeInp) = mean(maskedMov(icone,idx));
end

hold on;

fit_res = fit(iconeInp_log',spkRate','normcdf(x,mu,sigma)+b');
sig_log=[sig_log;fit_res.sigma];

if(~plot_fits)
plot(iconeInp_log,spkRate,col(jcone));
hold on;
else
plot(iconeInp_log',feval(fit_res,iconeInp_log'),col(jcone)) 
hold on;
end

end

% 
% 
% for jcone=[1:icone-1 ,icone+1:size(maskedMov,1)]
% %     pause
% thr = prctile(maskedMov(jcone,:),10);
% 
% for ibin=1:nbins+1
% X(ibin) = prctile(maskedMov(icone, maskedMov(jcone,:)>=thr),100*(ibin-1)/nbins);
% end
% 
% 
% iconeInp_log=[];
% spkRate=[];
% for iconeInp=1:nbins-1
%     
%     idx = (maskedMov(icone,:) > X(iconeInp) & maskedMov(icone,:) <= X(iconeInp+1) & maskedMov(jcone,:)>=thr);
% spkRate(iconeInp) = sum(binnedResponsesbigd(idx))/sum(idx);
% iconeInp_log(iconeInp) = mean(maskedMov(icone,idx));
% end
% 
% 
% 
% fit_res = fit(iconeInp_log',spkRate','normcdf(x,mu,sigma)+b');
% 
% if(~plot_fits)
% plot(iconeInp_log,spkRate,col(jcone));
% hold on;
% else
% plot(iconeInp_log',feval(fit_res,iconeInp_log'),col(jcone)) 
% hold on;
% end
% end
 