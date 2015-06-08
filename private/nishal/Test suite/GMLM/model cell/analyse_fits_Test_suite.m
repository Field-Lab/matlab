
kkx = fitGMLM.data_act.kx;
corr =zeros(nSU,nSU);

for isu=1:nSU
    for jsu=1:4
    r=corrcoef(kkx{isu}',cell_resp(:,jsu));
    corr(isu,jsu) = r(2,1);
    end
end

figure;
subplot(2,2,1)
scatter(kkx{1},kkx{2},4);
title('Raw dist of s.u. act');

subplot(2,2,2);
scatter(kkx{1}(binnedResponses~=0),kkx{2}(binnedResponses~=0),4);
title('spk trig distribution of s.u.');
subplot(2,2,3);

scatter(f(cell_resp(binnedResponses~=0,2)),f(cell_resp(binnedResponses~=0,3)),4)
title('True SU activation');

figure;
iidx=[1:1000];
plot(kkx{1}(iidx)/20);
hold on;
true = cell_resp(iidx,3) + cell_resp(iidx,4);
true = true.*(true>0);
plot(true);
legend('Model','True')


%see if it matches moment
pix_list = [1,2,5,4,3];
mov_len=size(maskedMov,2);
for ipix = pix_list
    for jpix=pix_list
    model_moment (ipix,jpix) = sum(maskedMov(ipix,:).*maskedMov(jpix,:).*fitGMLM.data_act.lam/120 )/ mov_len;
    rec_moment (ipix,jpix) = sum(maskedMov(ipix,binnedResponses~=0).*maskedMov(jpix,binnedResponses~=0))/mov_len;
    end
end
figure;plot(model_moment(:));hold on;plot(rec_moment(:));

% sta
model_sta = sum(maskedMov.*repmat(fitGMLM.data_act.lam/120,[size(maskedMov,1),1]),2)/mov_len;
rec_sta = sum(maskedMov(:,binnedResponses~=0),2)/mov_len;
