global server_path 
server_path = '/snle/lab/Experiments/Array/Analysis/';


%params.load_all = 1;
%TempDataRun = load_data('2008-03-25-4','rf-9-gf-nwpca-manual', params); 

if 1
    % initialize datarun
    datarun = load_data('/Analysis/gfield/2008-08-27-5/data005/data005');
    
    % load cell classification
    datarun = load_params(datarun,struct('verbose',1));
    
    % load STAs
    datarun = load_sta(datarun,struct('verbose',1));
    
    %temp_cell_types = order
    
    % set polarities
    params.cell_specs = {{6,9,11,13}, {7,8,10}}
    datarun = set_polarities(datarun, params);
    
    % compute summary frames
    datarun = save_sta_summaries(datarun,{6,7,8,9,10,11,13});
    
    % compute centers of mass (COM) for summary frames
    datarun = get_summary_com(datarun,{6,7,8,9,10,11,13});
    
    % get a quick look at ON parasol cells
    summary_portraits(datarun,{10},struct('plot_radius',5));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RGBtoLMS = [8.5667256e-7 3.2690937e-7 3.1051197e-8; 2.8230554e-6 3.051172e-6 1.5783603e-7; 6.9619523e-7 1.1124398e-6 1.8809839e-6];
RGBtoLMS = RGBtoLMS * 1e8;
LMStoRGB = inv(RGBtoLMS);
%LMStoRGB = [1796018.9 -187351.5 -13927.657; -1678713.5 513201.06 -15351.30; 1328066.3 -234171.1 545870.44];
%LMStoRGB = [1.0 -0.1043 -0.007755; -1.0 0.3057 -0.009145; 0.601 -0.429 1.0]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get data for LMS run
CellID = 1040 
datarunLMS = load_data('/Analysis/gfield/2008-08-27-5/data005/data005');
datarunLMS = load_params(datarunLMS,struct('verbose',1));

datarunLMS = load_sta(datarunLMS,struct('load_sta', [CellID]))'
CellNum = get_cell_numbers(datarunLMS, CellID)';

STASummary = one_sta_summary(datarunLMS.stas.stas{CellNum});

TempParams.tc_params.thresh = 20;
plot_time_course(datarunLMS.stas.stas{CellNum}, TempParams)

LMSTempImage = datarunLMS.stas.summaries{CellNum};
LMSTempImage = frame_image(LMSTempImage);
figure(2)
image(LMSTempImage)

LMSSTA = datarunLMS.stas.stas{CellNum};
figure(3)
hold on
plot(squeeze(LMSSTA(17,19,1,:)), 'r')
plot(squeeze(LMSSTA(17,19,2,:)), 'g')
plot(squeeze(LMSSTA(17,19,3,:)), 'b')
hold off


% convert to RGB prediction
LMSTempSTA = datarunLMS.stas.stas{CellNum};
[RowNum, ColNum, ColorNum,TimePnts] = size(LMSTempSTA);
TempSTA = LMSTempSTA;
PredictedRGBSTA = zeros(RowNum,ColNum,ColorNum,TimePnts);
Transform = RGBtoLMS';
for rw = 1:RowNum
    for cl = 1:ColNum
        for tpt = 1:TimePnts
            TempLMS = squeeze(TempSTA(rw,cl,:,tpt))';
            TempRGB = TempLMS * Transform;
            PredictedRGBSTA(rw,cl,:,tpt) = TempRGB';
            if (rw == 17) & (cl == 19)
                TempRGB
            end
        end
    end
end

TempFrame = squeeze(PredictedRGBSTA(:,:,:,29))
significant_stixels = threshold_sta(TempFrame, struct('thresh_type', 'peak','thresh', 0.2));
NumSigPix = length(find(significant_stixels > 0))
SigPix = 
for tpt = 1:TimePnts
    for clr = 1:NumColor
        tempclr = squeeze(PredictedRGBSTA(:,:,clr,tempclr)) .* significant_stixels;
        
        


% plot time course
TempParams.tc_params,
TempParams.tc_params.thresh = 5;
plot_time_course(PredictedRGBSTA)

figure(4)
hold on
plot(squeeze(PredictedRGBSTA(17,19,1,:)), 'r')
plot(squeeze(PredictedRGBSTA(17,19,2,:)), 'g')
plot(squeeze(PredictedRGBSTA(17,19,3,:)), 'b')
hold off

PredictedLMSSTASummary = one_sta_summary(PredictedLMSSTA);
% get summary
datarunLargeRGB = save_sta_summaries(datarunLargeRGB, [CellID]);
RGBTempImage = normalize_summary(RGBTempImage);
PredictedLMSScaled = frame_image(PredictedLMSSTASummary);
figure(6)
image(PredictedLMSScaled)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get data for RGB 10x10
CellID = 1040 
datarunLargeRGB = load_data('/Analysis/gfield/2008-08-27-5/data004/data004');
datarunLargeRGB = load_params(datarunLargeRGB,struct('verbose',1));

datarunLargeRGB = load_sta(datarunLargeRGB,struct('load_sta', [CellID]))'
CellNum = get_cell_numbers(datarunLargeRGB, CellID)';

STASummary = one_sta_summary(datarunLargeRGB.stas.stas{CellNum});
RGBTempImage = normalize_summary(STASummary);
ScaledRGBTempImage = frame_image(RGBTempImage);
figure(5)
image(ScaledRGBTempImage)

TempParams.tc_params.thresh = 15;
plot_time_course(datarunLargeRGB.stas.stas{CellNum}, TempParams)

TempSTA = datarunLargeRGB.stas.stas{CellNum};

figure(5)
hold on
plot(squeeze(TempSTA(17,19,1,:)), 'r')
plot(squeeze(TempSTA(17,19,2,:)), 'g')
plot(squeeze(TempSTA(17,19,3,:)), 'b')
hold off




%transform the entire STA
RGBTempSTA = datarunLargeRGB.stas.stas{CellNum};
[RowNum, ColNum, ColorNum,TimePnts] = size(RGBTempSTA);
TempSTA = RGBTempSTA;
PredictedLMSSTA = RGBTempSTA;
for rw = 1:RowNum
    for cl = 1:ColNum
        for tpt = 1:TimePnts
            TempRGB = squeeze(TempSTA(rw,cl,:,tpt))';
            TempLMS =  TempRGB * RGBtoLMS;
            PredictedLMSSTA(rw,cl,:,tpt) = TempLMS;
        end
    end
end
% plot time course
%figure(3)
TempParams.tc_params.thresh = 5;
plot_time_course(PredictedLMSSTA,TempParams)

PredictedLMSSTASummary = one_sta_summary(PredictedLMSSTA);
% get summary
datarunLargeRGB = save_sta_summaries(datarunLargeRGB, [CellID]);
RGBTempImage = normalize_summary(RGBTempImage);
PredictedLMSScaled = frame_image(PredictedLMSSTASummary);
figure(6)
image(PredictedLMSScaled)
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get data for RGB 1x1
datarunRGB = load_data('/Analysis/gfield/2008-08-27-5/data003/data003');
datarunRGB = load_params(datarunRGB,struct('verbose',1));

datarunRGB = load_sta(datarunRGB,struct('load_sta', [1038]))'
CellNum = get_cell_numbers(datarunRGB, 1038)';

datarunRGB = save_sta_summaries(datarunRGB, [1038]);

RGBTempImage = datarunRGB.stas.summaries{CellNum};
RGBTempImage = frame_image(RGBTempImage);
figure(2)
image(RGBTempImage)

BinnedSTA = rebin_sta(datarunRGB.stas.stas{CellNum}, 10);
figure(6)
hold on
plot(squeeze(BinnedSTA(17,19,1,:)), 'r')
plot(squeeze(BinnedSTA(17,19,2,:)), 'g')
plot(squeeze(BinnedSTA(17,19,3,:)), 'b')
hold off



BinnedSTA = rebin_sta(datarunRGB.stas.stas{CellNum}, 10);
size(BinnedSTA)

BinnedSTASummary = one_sta_summary(BinnedSTA);
ScaledBinnedSTASummary = frame_image(BinnedSTASummary);
figure(3)
image(ScaledBinnedSTASummary)



[RowNum, ColNum, ColorNum] = size(BinnedSTASummary);
TempSTA = BinnedSTASummary;
PredictedLMSSTASummary = BinnedSTASummary;
for rw = 1:RowNum
    for cl = 1:ColNum
        TempRGB = squeeze(TempSTA(rw,cl,:))';
        TempLMS = TempRGB * RGBtoLMS;
        PredictedLMSSTASummary(rw,cl,:) = TempLMS;
    end
end
PredictedLMSScaled = frame_image(PredictedLMSSTASummary);

figure(4)
image(PredictedLMSScaled)
        
        
        
        
        

