function fitGMLM_log = compute_op_NL_WN(WN_datafile,movie_xml,stim_length,cellID,user_STA_depth,ASM_link)


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

%% Find stimulus and response
% cellID=1531;

% user_STA_depth=30;
extract_movie_response2;

%% Load fits
data1 = load(ASM_link);
fitGMLM_log = data1.fitGMLM_log;

mask2=data1.mask;
mask=mask2(:);

%%
for nSU=1:length(fitGMLM_log)
fitGMLM = fitGMLM_log{nSU};

% compute NL
g = compute_op_NL(fitGMLM,maskedMovdd2,spksGen_hr);
fitGMLM_log{nSU}.NL_op = g;
end
end