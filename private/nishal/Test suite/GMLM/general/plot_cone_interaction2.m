function sig_log = plot_cone_interaction2(binnedResponsesbigd  , maskedMov ,ref_cone,cones1,cones2,plot_fits)

icone = ref_cone;sig_log=[];
col='rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr';
for jcone = cones1
col(jcone)='b';
end

for jcone = cones2
    col(jcone)='g';
end

%figure;
% subplot(3,1,1);

nbins=7;
for ibin=1:nbins+1
X(ibin) = prctile(maskedMov(icone,:),100*(ibin-1)/nbins);
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


    line_width=0.5;
    if(strcmp(col(jcone),'b')|strcmp(col(jcone),'g') )
    line_width=2;
    end
if(~plot_fits)

plot(iconeInp_log,spkRate,col(jcone),'LineWidth',line_width);
hold on;
else
fit_res = fit(iconeInp_log',spkRate','normcdf(x,mu,sigma)+b');
sig_log=[sig_log;fit_res.sigma];
plot(iconeInp_log',feval(fit_res,iconeInp_log'),col(jcone)) 
hold on;
end


end

spkRate=[];
for iconeInp=1:nbins-1
    idx = (maskedMov(icone,:) >= X(iconeInp) & maskedMov(icone,:) <= X(iconeInp+1));
spkRate(iconeInp) = sum(binnedResponsesbigd(idx))/sum(idx);
iconeInp_log(iconeInp) = mean(maskedMov(icone,idx));
end


if(~plot_fits)
plot(iconeInp_log,spkRate,'k','LineWidth',2);
hold on;
else
fit_res = fit(iconeInp_log',spkRate','normcdf(x,mu,sigma)+b');
sig_log=[sig_log;fit_res.sigma];
plot(iconeInp_log',feval(fit_res,iconeInp_log'),'k','LineWidth',2) 
hold on;
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
 