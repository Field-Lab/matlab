%% set up the simulation
clear all
%%
% stimulus parameters
stim.type = 'gaussian';
stim.len = [11520 1000];
stim.p = 1;

% weight and assignment parameters
sim.type = 'midget';
sim.C = 10;
[sim.I_sc datSim] = mkSubAssgn(sim.C,sim.type); % cone to subunit assignments
%sim = addNeighbors(dat,sim,0.5);
sim.S = size(sim.I_sc,1);
sim.C = size(sim.I_sc,2);
sim.J_s = ones(1,sim.S); % subunit to RGC assignments 

%sim.B_s = normMat(rand(1,sim.S),'L1',2)*5; % subunit to RGC weights
%sim.A_sc = normMat(ones(sim.S,sim.C).*sim.I_sc,'L1',2); % cone to subunit weights

load(fullfile('~/Dropbox/Projects/retinaSubunits/examples','midget-sim.mat'));
sim.B_s = fit.B_s;
sim.A_sc = fit.A_sc;
sim.f = fit.f;
sim.g = fit.g;
clear fit;

% nonlinearity parameters
%sim.f.type = 'loglinear';
%sim.f.p = [0 10 0.1 1];
%sim.f.type = 'loglinear';
%sim.f.p = [-2 2 -0.33 -0.03];
%sim.g.type = 'logexp1';
%sim.g.p = [0 1];

% generate stimuli
stim.C = sim.C;
[train test] = getStim(stim); % stimulus matrix

% cone non-linearity
%coneNonLin.type = 'logexp1';
%coneNonLin.p = [-2 1];
%inputTrain = (evalNonLin(train.X_ct,coneNonLin)-2);
%inputTest = (evalNonLin(test.X_ct,coneNonLin)-2);

% get mixing matrix
%sim.mix.type = 'e-neighborhood';
%sim.mix.p = [4 0.5];
%M = getMixMat(datSim.locs_c,sim.mix);
%M = eye(sim.C);

inputTrain = train.X_ct;
inputTest = test.X_ct;

%inputTrain = M*inputTrain;
%inputTest = M*inputTest;

% get the training responses
Y_st = evalNonLin((sim.I_sc.*sim.A_sc)*inputTrain,sim.f);
Z_t = evalNonLin((sim.J_s.*sim.B_s)*Y_st,sim.g);
train.R_t = poissrnd(Z_t);

% get the test responses
Y_st = evalNonLin((sim.I_sc.*sim.A_sc)*inputTest,sim.f);
Z_t = evalNonLin((sim.J_s.*sim.B_s)*Y_st,sim.g);
test.R_t = poissrnd(Z_t);

%% test greedy

datSim.nConesFit = sim.C;
datSim.polarity = 1;
fitGreedy(datSim,train,test,1);

%%
% fitting options
clear fit
fit.error = 'loglik';
fit.constraints.B_s = 'none';
fit.constraints.A_sc = 'L1';
fit.S = size(sim.I_sc,1);
fit.C = size(sim.I_sc,2);
%% try to estimate B_s assuming everything else

% stuff we'll assume
fit.I_sc = sim.I_sc;
fit.A_sc = sim.A_sc;
fit.J_s = sim.J_s;
fit.f = sim.f;

% stuff we are initializing
fit.B_s = normMat(rand(1,fit.S),'L1',2);
fit.g.type = 'logexp1';
fit.g.p = [0 1];
fit = fitB(train,fit,1);

%% try to estimate A_sc assuming everything else

% stuff we'll assume
fit.I_sc = sim.I_sc;
fit.J_s = sim.J_s;
fit.B_s = sim.B_s;
fit.g = sim.g;
fit.f = sim.f;

% stuff we are initializing
fit.A_sc = normMat(rand(sim.S,sim.C).*fit.I_sc,'L1',2);
fit = fitA(train,fit,1);

%% try to estimate g as a spline
% stuff we'll assume
fit.I_sc = sim.I_sc;
fit.A_sc = sim.A_sc;
fit.J_s = sim.J_s;
fit.J_s = sim.J_s;
fit.B_s = sim.B_s;
fit.g = sim.g;

% stuff we are initializing
fit.f.type = 'spline';
fit.f.init = 'linear';
fit.f.nknots = 8;
fit.f.smoothness = 3;
fit.f.extrap = 1;
fit.f.w = [];
fit = fitF(train,fit,0);
%% try to infer both
% stuff we'll assume
fit.I_sc = sim.I_sc;
%fit.I_sc = eye(fit.C);
fit.J_s = ones(1,size(fit.I_sc,1));

fit.S = size(fit.I_sc,1);
fit.C = size(fit.I_sc,2);

% stuff we are initializing
fit.f.type = 'linear';
fit.f.p = [0 1];
fit.g.type = 'logexp1';
fit.g.p = [0 1];
fit.B_s = normMat(ones(1,fit.S),'L1',2);
%fit.A_sc = normMat(ones(fit.S,fit.C).*fit.I_sc,'L1',2);
fit.A_sc = sim.A_sc;
fit.f = initSpline(fit.f);
fit = fitF(train,fit,1);
%fit.f.type = 'loglinear';
%fit.f.p = [-2 2 -0.33 -0.03];
%fit.f.p = [-1.5 0.5 5 -0.05];
%% try greedy

%%
plotModel(fit,datSim); drawnow;
nIter = 3;
for iIter = 1:nIter
    fit = fitB(train,fit,0);
    %fit = fitA(train,fit,0);
    fit = fitF(train,fit);
    plotModel(fit,datSim); drawnow;
end
%%
compareFit(dat,fit,sim)
%%
fit_tmp.g.type = 'exp';
fit_tmp.g.p = [0 1];
fit_tmp.f.type = 'rectLinear';
fit_tmp.f.p = [0 1];
%%
v = getPairwiseIntAll(train,test,fit_tmp);
%%
ind1 = 3; ind2 = 4;
out = getPairwiseInt(train,test,ind1,ind2,fit,0)
plotPairwiseInt(train,ind1,ind2,'gaussian',out,fit);
