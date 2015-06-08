function  respProb = getRespProb(fpath,patternNo,neuronId,stimAmp)
fname = ['elecResp_n' num2str(neuronId) '_p' num2str(patternNo) '.mat']; 
% Load elecRespFile
temp = load([fpath fname]);
elecResp = temp.elecResp; clear temp;
              
 plotResponseCurves = 0;                
[~, completeFit, erfErr] = fitToErf(elecResp,plotResponseCurves);  %#ok<ASGLU>
if stimAmp > max(completeFit(1,:))
    respProb = completeFit(2,end); 
    disp('Warning - response probability set by the maximum stimulation amplitude, not the one requested'); 
else
    respProb = completeFit(2,find(round(1000*completeFit(1,:)) == round(1000*stimAmp))); %#ok<FNDSB>
end
if isempty(respProb)
    disp('check function getRespProb')
    keyboard;
end
end