%This script generates a csv with success rates for each amplitude, in order to plot and analyze activation curve.

%File Names
%fn = 'E:\2015-12-18-0\data002\elecResp_n783_p50_r416_win.mat';
%outfn = 'E:\2015-12-18-0\data002\elecResp_n783_p50_r416.csv'
fn = '/Volumes/Analysis/2016-01-05-6/data003/elecResp_n2432_p181_r163.mat';
outfn = '/Volumes/Analysis/2016-01-05-6/data003/elecResp_n2432_p181_r163.csv';

%Load elecResp and get data
load(fn);
sa = elecResp.stimInfo.stimAmps;
sr = elecResp.analysis.successRates;
outm = [sa sr];


%write file
csvwrite(outfn, outm)
disp('Done!')
