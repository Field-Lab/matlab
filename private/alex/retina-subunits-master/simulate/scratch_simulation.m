load fit-2-subunits.mat
fit_2 = fit;

load fit-4-subunits.mat
fit_4 = fit;

load fit-8-subunits.mat
fit_8 = fit;

%% set up the simulation
fitOrig = fit_8;

% stimulus
stim.type = 'gaussian';
stim.grating = 'modulated';
stim.phase = pi/2;
stim.freq = myLogspace(log10(0.04),log10(0.1),10,10);
stim.len = [10000 5000];
stim.p = 1;
stim.C = 49;
[train test] = getStim(stim,dat); % stimulus matrix
sim = fitOrig;
%sim.g.p = 0;
%sim.f.p = 0;

% get the training responses
Y_st = evalNonLin((sim.I_sc.*sim.A_sc)*train.X_ct,sim.f);
Z_t = evalNonLin((sim.J_s.*sim.B_s)*Y_st,sim.g);
train.R_t = poissrnd(Z_t);

% get the test responses
Y_st = evalNonLin((sim.I_sc.*sim.A_sc)*test.X_ct,sim.f);
Z_t = evalNonLin((sim.J_s.*sim.B_s)*Y_st,sim.g);
test.R_t = poissrnd(Z_t);

%plot(Z_t)
%%

out(1) = evalFit(test,fit_2);
out(2) = evalFit(test,fit_4);
out(3) = evalFit(test,fit_8);
out(4) = evalFit(test,fit_16);
k = [2 4 8 16];

for i=1:4
nll(i) = out(i).negloglik;
r2(i) = out(i).r2;
end

base = min([nll]);
nll = nll-base;

subplot(2,1,1);
hold on
plot(k,r2,'k.-');
%drawHorzLine(r2_LN);
xlim([1 17]);
ylim([0 1]);

subplot(2,1,2);
hold on
plot(k,nll,'k.-');
%drawHorzLine(nll_LN);
xlim([1 17]);
ylim([0 20000]);

%%

load singleton-fit-2011-07-05-2-data002-3203.mat
fit_singleton = fit; clear fit;

for i=1:15; lik_sub(i) = evalFit(test,fits_sub{i}); end
lik_LN = evalFit(test,fit_LN);
lik_singleton = evalFit(test,fit_singleton);
hold on
plot(lik_sub,'k');
drawHorzLine(lik_LN,'b');
%drawHorzLine(lik_singleton,'r');
lik_sub_all(ifreq,iph,:) = lik_sub;
lik_LN_all(ifreq,iph,:) = lik_LN;
hold off

%%
clf
colormap(gray);
for i=1:size(test.X_ct,2);
    scatter(dat.locationsFit(:,1),dat.locationsFit(:,2),50,test.X_ct(:,i),'filled');
    set(gca,'CLim',[-2 2]);
    drawnow;
end
