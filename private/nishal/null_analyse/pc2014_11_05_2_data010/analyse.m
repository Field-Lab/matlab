addpath(genpath('../null_analyse/'));
startup_null_analyse_tenessee
%% Load cells and STAs
global_vars
datafile = '2014-11-05-2/data009';
type_name= cell(1,1);
type_name{1}=cell_params.type_name_inp;

datarun=load_data(datafile)
datarun=load_params(datarun)



cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

neuronPairsRefVsNew = crossIdentifyNeuronIDs('/Volumes/Analysis/2014-11-05-2/data009', '/Volumes/Analysis/2014-11-05-2/data010',InterestingCell_vis_id);
ref__cells=neuronPairsRefVsNew(:,2);


%%
