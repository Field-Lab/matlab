% Display the results for the case where the same electrode is used for
% stimulation and recording

% Load workspace
% Processed data for the same stimulating and recording electrode
% load(['/Users/grosberg/matlab/dataset_specific/GonzaloJune92015/'...
%     'sameRecandStimElec_29cases_10jun2015.mat']);


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
% load('/Users/grosberg/matlab/dataset_specific/workspace-axonal_activation_11jun2015.mat')
%%
for p = 541:555 % 555
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