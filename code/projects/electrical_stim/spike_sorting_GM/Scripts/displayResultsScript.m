% Display the results for the case where the same electrode is used for
% stimulation and recording

% Load workspace
% Processed data for the same stimulating and recording electrode
% load('/Volumes/Lab/Projects/electrical_stim/GM-sorting-validation/GonzaloJune92015/sameRecandStimElec_29cases_10jun2015.mat')

%%
for p = 1:size(NeuronList,1)
    Gibbs = sol(p).Gibbs;
    GibbsNoDelete = sol(p).GibbsNoDelete;
    initial = sol(p).initial ;
    input = sol(p).input;
    Log = sol(p).Log;
    %times(p)=toc;
    time=times(p);
    %      DisplayResults(input,GibbsNoDelete,Log,1)
    displayResults_compare(input, Gibbs, Log, 1,1)
end
 
%%
% load('/Volumes/Lab/Projects/electrical_stim/GM-sorting-validation/2015-06-algorithm-results/workspace-axonal_activation_11jun2015.mat')
%% Load the cases where a different 
for p = 450; %541:555 % p can range from 450 to 555
    Gibbs = solAxonB(p).Gibbs;
    GibbsNoDelete = solAxonB(p).GibbsNoDelete;
    initial = solAxonB(p).initial ;
    input = solAxonB(p).input;
    Log = solAxonB(p).Log;
   
    %      DisplayResults(input,GibbsNoDelete,Log,1)
    displayResults_compare(input, Gibbs, Log, 1)
end

plotResponseCurves = 0 ; 
[threshold, completeFit, erfErr] = fitToErf(elecResp,plotResponseCurves); 
%%
for p = 1:20 % 555
    Gibbs = solAxon(p).Gibbs;
    GibbsNoDelete = solAxon(p).GibbsNoDelete;
    initial = solAxon(p).initial ;
    input = solAxon(p).input;
    Log = solAxon(p).Log;
    displayResults_compare(input, Gibbs, Log, 1)
end

%% For axonal activation solutions
values = [0;0;0;0;0;0;0;0]; 
for p = [1:200] % 555
    Gibbs = solAxon(p).Gibbs;
    GibbsNoDelete = solAxon(p).GibbsNoDelete;
    initial = solAxon(p).initial ;
    input = solAxon(p).input;
    Log = solAxon(p).Log;
    values = algHumCmp(input,Gibbs,1,values);
    disp(num2str(p)); 
end

percentFalsePositive = values(4)/values(6);
percentTrueNegative = values(5)/values(6);
percentTruePositive = values(1)/values(3);
percentFalseNegatives = values(2)/values(3);
percentAgreement = values(7)/values(8);
totalTrials = values(8); 
%%
groupnames = {'overall success','false positive','true negative',...
    'true positive','false negative'}; 
bw_title = 'success rates'; 
bw_xlabel ='type'; 
bw_ylabel = 'percentage';
gridstatus = 'xy';
bw_legend = {'red','blue'};
my_errors = [sqrt(totalTrials) sqrt(totalHumanNegatives) ...
    sqrt(totalHumanNegatives) sqrt(totalHumanPositives) sqrt(totalHumanPositives)]; 

all_values = [percentAgreement; percentFalsePositive; percentTrueNegative; ...
    percentTruePositive; percentFalseNegatives]; 
bh = barweb([all_values],[my_errors;]', [], groupnames, bw_title, ...
    bw_xlabel, bw_ylabel, [], gridstatus, bw_legend);

 figure; bar(all_values); 
 ylabel('percent accuracy');