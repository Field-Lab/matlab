%% find_data2.m
% As pattern number, neuron numebr etc provided 

% 
% stim 4
% cell_no = 197;
% patternNo = 27;

% stim5 
%cell_no= 436
%patternNo = 76
% istim=1;

%%
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile('/Volumes/Analysis/2012-09-24-3/data007/data007.ei');
ids = eiFile.getIDList();

ei = eiFile.getImage(cell_no,0);  % (id, errorType (Standard Deviation of the Mean: 0, Variance of the Mean: 1))
% for standard version of matlab, no second argument (standard deviation vs. variance not an option)

% returns ei(average:1 error:2, electrode + 1, time index)
% 3 dimensional array with first dimension referring to whether value is the average voltage (1) or the error (2), second dimension referring to electrode number, and % third dimension referring to point in time

maxElectrode = eiFile.getMaxElectrode(ei);

%%
pathToAnalysisData = '/Volumes/Analysis/nishal/data008/';

movieNos = [];
patternNoString = ['p' num2str(patternNo)];
files = dir([pathToAnalysisData patternNoString]);

for i = 1:length(files)
   if strfind(files(i).name, patternNoString) == 1
       mIndices = strfind(files(i).name, 'm');
       movieNos = [movieNos str2double(files(i).name(mIndices(end)+1:end))]; %#ok<AGROW>
   end
end
movieNos = sort(movieNos);

% %%
% lam =2.5;%1
% eps=0.00001;%0.05
% gam=1;%1
% mu=100000000;%100000
% movieNo = 306;%movieNos(2);
% idxs=[1:length(movieNos)];
% movie_idx=idxs(movieNos==movieNo);
% fcn_real_anal(istim,movie_idx,movieNo, patternNo,cell_no, ei,maxElectrode,mu,lam,eps,gam);


%%
for mu=100000000%[100]
for lam=2.5%[0.1,0.3,0.5,0.7,1,1.2,1.5]
for eps=0.00001%[1,0.5,0.1,0.05,0.01]
for gam=1%0.5:0.5:1.5
parfor imov=1:length(movieNos)
    try
movieNo = movieNos(imov);
fcn_real_anal(istim,imov,movieNo, patternNo,cell_no, ei,maxElectrode,mu,lam,eps,gam);
    catch
    end
 end
end
end
end
end

%%

load(sprintf('/Volumes/Analysis/2012-09-24-3/data008/elecResp_n%d_p%d.mat',cell_no,patternNo));

for threshold_resp=0.1:0.1:0.8
prob_resp_log=0*[1:length(movieNos)];


for imov=1:length(movieNos)
load(sprintf('/Volumes/Analysis/nishal/2012-09-24-3/stim%d/mov%d.mat',istim,movieNos(imov)));

indx=[1:length(data_event)];
dat_pts = indx(data_event==1);
analysis_window = 0*x;
for i=1:100
    analysis_window(dat_pts+i-1)=1;
end
dat_windowed = y.*x.*analysis_window;

% for each stimulation
resp_cnt=0;
for iistim=1:sum(data_event)
    if(sum(dat_windowed(dat_pts(iistim):dat_pts(iistim)+100)>threshold_resp)>0)
    resp_cnt=resp_cnt+1;
    end
end
prob_resp=resp_cnt/sum(data_event);

prob_resp_log(imov) = prob_resp;


% Correct probability distribution
%prob_resp_cnt_human = sum(elecStim.latencies{stimulation_data}>0)/length(elecStim.latencies{stimulation_data});
%prob_resp_log_human(stimulation_data) = prob_resp_cnt_human;
end
fi = figure;
plot([1:length(movieNos)],prob_resp_log,'b');
hold on;
plot([1:length(movieNos)],elecResp.analysis.successRates,'r');
ylim([0,1.5])
%plot([1:69],prob_resp_log_human,'r');
legend('CBP','human','Location','NorthWest');
title(sprintf('Threshold %f',threshold_resp));
print(fi,sprintf('/Volumes/Analysis/nishal/2012-09-24-3/stim%d/prob_resp_curve_th_%f.eps',istim,threshold_resp),'-depsc');
end



%% 