function  LL = getLL_const_var_2015_03_09_2(spkCondColl,inputType)

for icond = 1:length(spkCondColl)
    if(inputType==1)
spks = makeSpikeMat(spkCondColl(icond).spksColl,1/120,1272);
    else
spks = spkCondColl(icond).spksColl;    
    end
    
binlen = 12;
%[t,rt] = psth_calc(spks,binlen,'overlap'); %mean(spks,1);
%rt = [zeros(1,binlen/2),rt,zeros(1,-1+(binlen/2))];
rt = mean(spks,1);
rt(rt==0) = 0.000001;
r = mean(rt)*ones(size(rt));

tms=(rt>0);
rt =rt(tms);
r=r(tms);
spks=spks(:,tms);

LL_s = 0;
LL_b = 0 ;
for itrial =1:size(spks,1)
LL_s = LL_s + -sum(r) + spks(itrial,:)*log(r)' ; % - sum(log(factorial(spks(itrial,:))));
LL_b = LL_b + -sum(rt) + spks(itrial,:)*log(rt)'; % - sum(log(factorial(spks(itrial,:))));
end
LL_s = LL_s / sum(spks(:));
LL_b = LL_b / sum(spks(:));

LL(icond,1) = LL_s; LL(icond,2)=LL_b;
end


end

