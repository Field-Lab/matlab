function [nlike_full,nlike_diff,diffll]=compare_lam_likelihood(lam1,lam2,binnedResponses_global)

tsp = find(binnedResponses_global~=0);
y_tsp = binnedResponses_global(binnedResponses_global~=0)';
dt= 1/1200;

lam = lam1;
likelihood = sum(-lam*dt) + sum(y_tsp.*log(lam(tsp)));
fval1=-likelihood/length(binnedResponses_global);

lam = lam2;
likelihood = sum(-lam*dt) + sum(y_tsp.*log(lam(tsp)));
fval2=-likelihood/length(binnedResponses_global);
nlike_full = [fval1,fval2];


%% difference
lam_diff = abs(lam1-lam2);
cutoff = prctile(lam_diff,80);
tcutoff = lam_diff>cutoff;
ttsp = binnedResponses_global~=0;

tsp_cutoff = ttsp & tcutoff';
y_tsp_cutoff = binnedResponses_global(tsp_cutoff)';

lam = lam1;
likelihood = sum(-lam(tcutoff)*dt) + sum(y_tsp_cutoff.*log(lam(tsp_cutoff)));
fval_cutoff1=-likelihood/sum(tcutoff);

lam = lam2;
likelihood = sum(-lam(tcutoff)*dt) + sum(y_tsp_cutoff.*log(lam(tsp_cutoff)));
fval_cutoff2=-likelihood/sum(tcutoff);

nlike_diff=[fval_cutoff1,fval_cutoff2];


%% difference LL  


l1 = lam1(tcutoff); l2 =lam2(tcutoff);
y =  binnedResponses_global(tcutoff);
 ll1 = -l1*dt + y'.*log(l1);
 ll2 = -l2*dt + y'.*log(l2);
 
 diffll = ll1-ll2;

end