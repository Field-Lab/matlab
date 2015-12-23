
load('/Analysis/lauren/2012-01-27-3/data002/elecResp_n93_p5_w50.mat')

latencies = elecResp.analysis.latencies{23};

latenciesForHist = latencies(latencies~=0)/20;

binEdges = (0:40)/40;
binCounts = histc(latenciesForHist, binEdges);
binCounts = binCounts(1:end-1); %last bin only includes latencies == tempMinEndMs;

data = zeros(2,length(binEdges)-1);
data(1,:) = binEdges(1:end-1);
data(2,:) = binCounts;

t0Guess = 0.2;
tauGuess = 0.2;
n = 2;

[estParams alpha] = impRespFitter(data, t0Guess, tauGuess, n, 'makePlot', true);
%[estParams alpha] = impRespFitter(data, tauGuess, n, 'makePlot', true);



%% messing with the parameters

% n = 3;
% t0 = 2;
% tau = 10;

t0 = estParams(1);
tau = estParams(2);

t = 0:0.001:2;

t_sc = (t - t0)/tau; %scaled/shifted time values
p1 = exp(-n*t_sc);
p2 = t_sc.^n;
p = p1.*p2;


figure
hold on
plot(t,alpha*p, 'k-')
axis manual
plot(t,alpha*p1,'r-')
plot(t,alpha*p2, 'b-')


% t = -0.03:0.01:2;
% p_new = (t.*exp(-t/tau)).^n;
% 
% plot(t,p_new,'k-')