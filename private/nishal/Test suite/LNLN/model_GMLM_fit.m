

 mask2 = logical(ones(size(movie,1),size(movie,2)));
 maskedMov= filterMov(movie,mask2,squeeze(model.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];

 binnedResponses = response';

 filteredStimDim =size(maskedMov,1);
 
 %  EM like Max Expected Likelihood .. 
 interval=1;
 %[fitGMLM,output] = fitGMLM_MEL_EM(binnedResponses,maskedMov2,8,4,interval);   
 
 fitGMLM1=cell(15,1); 
 for nSU=1:15
 nSU
 fitGMLM_log = cell(50,1);
 fval_log = zeros(50,1);
 for ifit = 1:1
     ifit
     close all
 [fitGMLM,f_val] = fitGMLM_EM_bias(binnedResponses,maskedMov,filteredStimDim,nSU,interval); 
 
 fitGMLM_log{ifit} = fitGMLM;
 fval_log(ifit) = f_val;
 
 end
 fitGMLM1{nSU} = fitGMLM_log;
 %save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/stix 16, 90 min/fit_nSU_%d.mat',nSU),'fitGMLM_log','fval_log','model','mask2');
 end
