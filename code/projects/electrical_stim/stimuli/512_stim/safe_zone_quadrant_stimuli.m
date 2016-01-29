% Script used in 2015-10-29 safe zone tests.  
% Designed to measure the combination of electrode pairs that are
% perpendicular or parallel to the axon bundle to determine a model for
% safe zone predictions. 
% SET TRIGGER INTERVAL TO 0.5 SECONDS

clear all;

% as seen in vision-interactive (EI plot):
% 
% 2   |   1
%-----------
% 3   |   4

%% Define inputs
include_single_elecs = true;
both_polarity_combs = false;
quadrant = 1234; % 1, 2, 3, 4, 12, or 34 (12 does something in the middle of 1&2, 34 does 3 and 4 together)
elec_spacing = 60;
same_polarity = false; % Set to true to use only negative pairs. 
use_ratios = false;
scale_factor = 1; %Set scale factor to determine 2-elec ratio to use. A positive value gives
pairOrientation = 'downleft'; %horizontal, downleft, downright, vertical (horizontal must be for 60 µm and vertical must be for 30 µm)
secondOrientation = 'horizontal'; 
delayInMs = 7.5; %interval between pulses
saveFiles = 0; %Set to 1 to save stimulus files, 0 for testing
saveName = 'safezone_quad1234_2016_01_05_4'; %Descriptive name for the stimulus files 

%%
if elec_spacing == 60
    elec_coords = electrode_positions(512);
elseif elec_spacing == 30
    elec_coords = electrode_positions(519);
else
    error('invalid electrode spacing value')
end

if elec_spacing == 60
    if quadrant == 1
        electrodes = find(elec_coords(:,1)>=0 & elec_coords(:,2)>=0);
    elseif quadrant == 2;
        electrodes = find(elec_coords(:,1)<0  & elec_coords(:,2)>=0);
    elseif quadrant == 3;
        electrodes = find(elec_coords(:,1)<0 & elec_coords(:,2)<0);
    elseif quadrant == 3.1
        electrodes = find(elec_coords(:,1) <40 & elec_coords(:,1)>-940 ...
            & elec_coords(:,2)<0); 
    elseif quadrant == 34;
        electrodes = find(elec_coords(:,1)<472.5 & elec_coords(:,1)>-472.5...
            & elec_coords(:,2)<0);
    elseif quadrant == 12;
        electrodes = find(elec_coords(:,1)<472.5 & elec_coords(:,1)>-472.5...
            & elec_coords(:,2)>0);
    elseif quadrant == 1234;
        electrodes = find(elec_coords(:,1)<472.5 & elec_coords(:,1)>-472.5...
            & elec_coords(:,2)>-225 & elec_coords(:,2)<225);   
    else
        electrodes = find(elec_coords(:,1)>=0 & elec_coords(:,2)<0);
    end
elseif elec_spacing == 30 %leave out electrodes at x or y ~= 0
    if quadrant == 1
        electrodes = find(elec_coords(:,1)>0.1  & elec_coords(:,2)>0.1);
    elseif quadrant == 2;
        electrodes = find(elec_coords(:,1)<-0.1 & elec_coords(:,2)>0.1);
    elseif quadrant == 3;
        electrodes = find(elec_coords(:,1)<-0.1 & elec_coords(:,2)<-0.1);
    else
        electrodes = find(elec_coords(:,1)>0.1  & elec_coords(:,2)<-0.1);
    end
end

%%

nSamples  =10000; %0.5 s
timeShift =0;
delay     = delayInMs*20;

pPerChunk = floor(nSamples/delay);

array  = zeros(length(electrodes), 0);
elec_1 = zeros(1,0); %keeps track of electrode 1 (for ordering)
stim_type = zeros(1,0);

for ii = 1:length(electrodes)
    rel_coords = zeros(size(elec_coords));
    rel_coords(:,1) = elec_coords(:,1) - elec_coords(electrodes(ii),1);
    rel_coords(:,2) = elec_coords(:,2) - elec_coords(electrodes(ii),2);

    switch lower(pairOrientation)
        case{'horizontal'}
            if elec_spacing == 60
                elec_2 = find(rel_coords(:,1)>55 & rel_coords(:,1)<65 & rel_coords(:,2) == 0);
            else
                error('horizontal orientation only allowed for 60 micron spacing')
            end
        case{'vertical'}
            if elec_spacing == 30
                elec_2 = find(abs(rel_coords(:,1))<0.5 & rel_coords(:,2)>28 & rel_coords(:,2)<32);
            else
                error('vertical orientation only allowed for 30 micron spacing')
            end
        case{'downright'}
            if elec_spacing == 60
                elec_2 = find(rel_coords(:,1)>25 & rel_coords(:,1)<35 & rel_coords(:,2)>-65 & rel_coords(:,2)<-55);
            else
                elec_2 = find(rel_coords(:,1)>-32 & rel_coords(:,1)<-28 & rel_coords(:,2)>13 & rel_coords(:,2)<17);
            end
        case{'downleft'}
            if elec_spacing == 60
                elec_2 = find(rel_coords(:,1)>25 & rel_coords(:,1)<35 & rel_coords(:,2)>55 & rel_coords(:,2)<65);
            else
                elec_2 = find(rel_coords(:,1)>28 & rel_coords(:,1)<32 & rel_coords(:,2)>13 & rel_coords(:,2)<17);
            end
        otherwise
            error('invalid electrode pair orientation')
    end % Switch statement testing pair orientation
    if isempty(elec_2)
        disp(['excluding electrode ' num2str(electrodes(ii)) ' because there is no electrode in chosen direction'])
    elseif ~ismember(elec_2, electrodes)
        disp(['excluding electrode ' num2str(electrodes(ii)) ' because other electrode isn''t in this quadrant'])
    else
        array = [array zeros(length(electrodes),1)];
        array(ii, end) = 1;
        array(electrodes==elec_2, end) = -1;
        elec_1 = [elec_1 electrodes(ii)];
        stim_type = [stim_type 1];

        if both_polarity_combs
            array = [array zeros(length(electrodes),1)];
            array(ii, end) = -1;
            array(electrodes==elec_2, end) = 1;
            elec_1 = [elec_1 electrodes(ii)];
            stim_type = [stim_type 2];
        end
    end
end

if include_single_elecs
    array = [array eye(length(electrodes))];
    elec_1 = [elec_1 electrodes'];
    stim_type = [stim_type 3*ones(1,length(electrodes))];
end

%% order stimuli so that subsequent elec_1s are at least 3 electrode
% spacings apart

% break into 8 groups and select randomly from a sequence of groups
if elec_spacing == 60
    x_coords_range = [min(elec_coords(electrodes,1)) max(elec_coords(electrodes,1))];
    boundaries = x_coords_range(1)+(x_coords_range(2)-x_coords_range(1))*[0.125 0.25 0.375 0.5 0.625 0.75 0.875];
    
    e_groups{1} = electrodes(elec_coords(electrodes,1)<boundaries(1));
    e_groups{8} = electrodes(elec_coords(electrodes,1)>=boundaries(7));
    for ii = 2:7
        e_groups{ii} = electrodes(elec_coords(electrodes,1)>=boundaries(ii-1)...
            & elec_coords(electrodes,1)<boundaries(ii)); %#ok<SAGROW>
    end
else
   [~, sortInd] = sort(elec_coords(electrodes,1));
   for ii = 1:8
       e_groups{ii} = electrodes(sortInd(15*(ii-1)+1:15*ii));
   end
end

%%

if elec_spacing == 60
%     if quadrant == 34
%         pattern_order1 = zeros(1,256);
%         pattern_order2 = zeros(1,256);
%     else
        pattern_order1 = zeros(1,128);
        pattern_order2 = zeros(1,128);
%     end
else
    pattern_order1 = zeros(1,120);
    pattern_order2 = zeros(1,120);
end
    
e_group_order = [1 5 2 6 3 7 4 8];
for ii = 1:8
    iGroup = e_group_order(ii);
    pattern_order1(ii:8:end) = e_groups{iGroup}(randperm(length(e_groups{iGroup})));
    pattern_order2(ii:8:end) = e_groups{iGroup}(randperm(length(e_groups{iGroup})));
end


%%
% first two chunks: single bipolar polarity combination

if elec_spacing == 60
    for ii = 1:32
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{1}(ii) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{1}(ii+32) = pInd;
        end
    end
    
    for ii = 33:64
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{2}(ii-32) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{2}(ii) = pInd;
        end
    end
    
    for ii = 65:96
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{3}(ii-64) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{3}(ii-32) = pInd;
        end
    end
    
    for ii = 97:128
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{4}(ii-96) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{4}(ii-64) = pInd;
        end
    end
    
    if both_polarity_combs

        for ii = 1:32
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{5}(ii) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{5}(ii+32) = pInd;
            end
        end
        
        for ii = 33:64
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{6}(ii-32) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{6}(ii) = pInd;
            end
        end
        
        for ii = 65:96
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{7}(ii-64) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{7}(ii-32) = pInd;
            end
        end
        for ii = 97:128
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{8}(ii-96) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{8}(ii-64) = pInd;
            end
        end
    end
    
    if include_single_elecs
        
        for ii = 1:32
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{5+4*both_polarity_combs}(ii) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{5+4*both_polarity_combs}(ii+32) = pInd;
            end
        end
        
        for ii = 33:64
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{6+4*both_polarity_combs}(ii-32) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{6+4*both_polarity_combs}(ii) = pInd;
            end
        end
        
        for ii = 65:96
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{7+4*both_polarity_combs}(ii-64) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{7+4*both_polarity_combs}(ii-32) = pInd;
            end
        end
        
        for ii = 97:128
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{8+4*both_polarity_combs}(ii-96) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{8+4*both_polarity_combs}(ii-64) = pInd;
            end
        end
    end
else
    for ii = 1:30
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{1}(ii) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{1}(ii+30) = pInd;
        end
    end
    
    for ii = 31:60
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{2}(ii-30) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{2}(ii) = pInd;
        end
    end
    
    for ii = 61:90
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{3}(ii-60) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{3}(ii-30) = pInd;
        end
    end
    
    for ii = 91:120
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{4}(ii-90) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{4}(ii-60) = pInd;
        end
    end
    
    if both_polarity_combs
        for ii = 1:30
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{5}(ii) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{5}(ii+30) = pInd;
            end
        end
        
        for ii = 31:60
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{6}(ii-30) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{6}(ii) = pInd;
            end
        end
        
        for ii = 61:90
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{7}(ii-60) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{7}(ii-30) = pInd;
            end
        end
        for ii = 91:120
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{8}(ii-90) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{8}(ii-60) = pInd;
            end
        end
    end
    
    if include_single_elecs
        
        for ii = 1:30
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{5+4*both_polarity_combs}(ii) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{5+4*both_polarity_combs}(ii+30) = pInd;
            end
        end
        
        for ii = 31:60
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{6+4*both_polarity_combs}(ii-30) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{6+4*both_polarity_combs}(ii) = pInd;
            end
        end
        
        for ii = 61:90
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{7+4*both_polarity_combs}(ii-60) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{7+4*both_polarity_combs}(ii-30) = pInd;
            end
        end
        
        for ii = 91:120
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{8+4*both_polarity_combs}(ii-90) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{8+4*both_polarity_combs}(ii-60) = pInd;
            end
        end
    end
end

%% check that each stimulus is applied twice
times_applied = zeros(1,size(array,2));
for ii = 1:length(pattern_order_all)
    for jj =1:length(pattern_order_all{ii})
        if pattern_order_all{ii}(jj)~=0
            times_applied(pattern_order_all{ii}(jj)) = times_applied(pattern_order_all{ii}(jj))+1;
        end
    end
end

if ~all(times_applied == 2)
    error('each pattern should be applied twice...so something''g wrong!')
end


%%
MovieChunks = [];
for ii = 1:length(pattern_order_all)
    thisPattern = pattern_order_all{ii};
    Times = (0:delay:delay*(length(thisPattern)-1)) + timeShift;
    
    %get rid of 0 patterns
    for jj = 1:length(thisPattern)
        Times(thisPattern == 0) = [];
        thisPattern(thisPattern == 0) = [];
    end
    
    if any(Times > nSamples)
        error('stimulus time exceeds movie chunk time')
    end
        
    Chunk=NS_MovieChunkGenerationForExperiment(Times,nSamples,thisPattern);
    MovieChunks = [MovieChunks Chunk];
end
MovieChunksFile=[length(pattern_order_all) MovieChunks]; %only one movie


if same_polarity
    figure; subplot(1,2,1); imagesc(array); title('opp polarity pairs')
    array(find(array == -1)) = 1; 
    subplot(1,2,2); imagesc(array); title('same polarity pairs');
end

if use_ratios
    figure; subplot(1,2,1); imagesc(array); title('opp polarity pairs')
    if scale_factor<0
        array(find(array == 1)) = -scale_factor;
    else
        array(find(array == -1)) = -scale_factor; 
    end
    subplot(1,2,2); imagesc(array); title('ratio');
end
%% Save files
% cd /Users/grosberg/Desktop/; 
if saveFiles
    fid = fopen([saveName '_el'],'wb','ieee-le.l64');
    fwrite(fid,electrodes,'int32');
    fclose(fid);
    
    fid = fopen([saveName '_mv'],'wb','ieee-le.l64');
    fwrite(fid,MovieChunksFile,'int32');
    fclose(fid);
end
figure; for mm = 1:size(pattern_order_all,2); hold on; plot(pattern_order_all{mm},'-x'); end
for mm = 1:size(pattern_order_all,2); hold on; plot(diff(pattern_order_all{mm}),'-x'); end

clear pattern_order1 pattern_order2 pattern_order_all rel_coords ;
clear thisPattern times_applied Times;
%

%% now do the horizontal approach. 
electrodes1 = electrodes; 
addon = zeros(length(electrodes),11); 
addon(1:11,1:11) = eye(11);
array1 = [array(:,1:2:end) addon]; 
clear pattern_order_all; 
clear times_applied; 
%% Define inputs
% include_single_elecs = false;
% both_polarity_combs = false;
% quadrant = 3; % 1, 2, 3, 4, or 34 (34 does 3 and 4 together)
% elec_spacing = 60;
same_polarity = false; % Set to true to use only negative pairs. 
% use_ratios = false;
% scale_factor = -0.5; %Set scale factor to determine 2-elec ratio to use. A positive value gives
pairOrientation = secondOrientation; %horizontal, downleft, downright, vertical (horizontal must be for 60 µm and vertical must be for 30 µm)


%%
if elec_spacing == 60
    elec_coords = electrode_positions(512);
elseif elec_spacing == 30
    elec_coords = electrode_positions(519);
else
    error('invalid electrode spacing value')
end


if elec_spacing == 60
    if quadrant == 1
        electrodes = find(elec_coords(:,1)>=0 & elec_coords(:,2)>=0);
    elseif quadrant == 2;
        electrodes = find(elec_coords(:,1)<0  & elec_coords(:,2)>=0);
    elseif quadrant == 3;
        electrodes = find(elec_coords(:,1)<0 & elec_coords(:,2)<0);
    elseif quadrant == 3.1
        electrodes = find(elec_coords(:,1) <40 & elec_coords(:,1)>-940 ...
            & elec_coords(:,2)<0); 
    elseif quadrant == 34;
        electrodes = find(elec_coords(:,1)<472.5 & elec_coords(:,1)>-472.5...
            & elec_coords(:,2)<0);
    elseif quadrant == 12;
        electrodes = find(elec_coords(:,1)<472.5 & elec_coords(:,1)>-472.5...
            & elec_coords(:,2)>0);
    elseif quadrant == 1234;
        electrodes = find(elec_coords(:,1)<472.5 & elec_coords(:,1)>-472.5...
            & elec_coords(:,2)>-225 & elec_coords(:,2)<225);
    else
        electrodes = find(elec_coords(:,1)>=0 & elec_coords(:,2)<0);
    end
elseif elec_spacing == 30 %leave out electrodes at x or y ~= 0
    if quadrant == 1
        electrodes = find(elec_coords(:,1)>0.1  & elec_coords(:,2)>0.1);
    elseif quadrant == 2;
        electrodes = find(elec_coords(:,1)<-0.1 & elec_coords(:,2)>0.1);
    elseif quadrant == 3;
        electrodes = find(elec_coords(:,1)<-0.1 & elec_coords(:,2)<-0.1);
    else
        electrodes = find(elec_coords(:,1)>0.1  & elec_coords(:,2)<-0.1);
    end
end

%%

nSamples=10000; %0.5 s

timeShift=0;
delay=delayInMs*20;

pPerChunk = floor(nSamples/delay);

array = zeros(length(electrodes), 0);
elec_1 = zeros(1,0); %keeps track of electrode 1 (for ordering)
stim_type = zeros(1,0);

for ii = 1:length(electrodes)
    rel_coords = zeros(size(elec_coords));
    rel_coords(:,1) = elec_coords(:,1) - elec_coords(electrodes(ii),1);
    rel_coords(:,2) = elec_coords(:,2) - elec_coords(electrodes(ii),2);

    switch lower(pairOrientation)
        case{'horizontal'}
            if elec_spacing == 60
                elec_2 = find(rel_coords(:,1)>55 & rel_coords(:,1)<65 & rel_coords(:,2) == 0);
            else
                error('horizontal orientation only allowed for 60 micron spacing')
            end
        case{'vertical'}
            if elec_spacing == 30
                elec_2 = find(abs(rel_coords(:,1))<0.5 & rel_coords(:,2)>28 & rel_coords(:,2)<32);
            else
                error('vertical orientation only allowed for 30 micron spacing')
            end
        case{'downright'}
            if elec_spacing == 60
                elec_2 = find(rel_coords(:,1)>25 & rel_coords(:,1)<35 & rel_coords(:,2)>-65 & rel_coords(:,2)<-55);
            else
                elec_2 = find(rel_coords(:,1)>-32 & rel_coords(:,1)<-28 & rel_coords(:,2)>13 & rel_coords(:,2)<17);
            end
        case{'downleft'}
            if elec_spacing == 60
                elec_2 = find(rel_coords(:,1)>25 & rel_coords(:,1)<35 & rel_coords(:,2)>55 & rel_coords(:,2)<65);
            else
                elec_2 = find(rel_coords(:,1)>28 & rel_coords(:,1)<32 & rel_coords(:,2)>13 & rel_coords(:,2)<17);
            end
        otherwise
            error('invalid electrode pair orientation')
    end
    if isempty(elec_2)
        disp(['excluding electrode ' num2str(electrodes(ii)) ' because there is no electrode in chosen direction'])
    elseif ~ismember(elec_2, electrodes)
        disp(['excluding electrode ' num2str(electrodes(ii)) ' because other electrode isn''t in this quadrant'])
    else
        array = [array zeros(length(electrodes),1)];
        array(ii, end) = 1;
        array(electrodes==elec_2, end) = -1;
        elec_1 = [elec_1 electrodes(ii)];
        stim_type = [stim_type 1];

        if both_polarity_combs
            array = [array zeros(length(electrodes),1)];
            array(ii, end) = -1;
            array(electrodes==elec_2, end) = 1;
            elec_1 = [elec_1 electrodes(ii)];
            stim_type = [stim_type 2];
        end
    end
end

if include_single_elecs
    array = [array eye(length(electrodes))];
    elec_1 = [elec_1 electrodes'];
    stim_type = [stim_type 3*ones(1,length(electrodes))];
end

%% order stimuli so that subsequent elec_1s are at least 3 electrode
% spacings apart

% break into 8 groups and select randomly from a sequence of groups
if elec_spacing == 60
    x_coords_range = [min(elec_coords(electrodes,1)) max(elec_coords(electrodes,1))];
    boundaries = x_coords_range(1)+(x_coords_range(2)-x_coords_range(1))*[0.125 0.25 0.375 0.5 0.625 0.75 0.875];
    
    e_groups{1} = electrodes(elec_coords(electrodes,1)<boundaries(1));
    e_groups{8} = electrodes(elec_coords(electrodes,1)>=boundaries(7));
    for ii = 2:7
        e_groups{ii} = electrodes(elec_coords(electrodes,1)>=boundaries(ii-1) & elec_coords(electrodes,1)<boundaries(ii)); %#ok<SAGROW>
    end
else
   [~, sortInd] = sort(elec_coords(electrodes,1));
   for ii = 1:8
       e_groups{ii} = electrodes(sortInd(15*(ii-1)+1:15*ii));
   end
end

%%

if elec_spacing == 60
%     if quadrant == 34
%         pattern_order1 = zeros(1,256);
%         pattern_order2 = zeros(1,256);
%     else
        pattern_order1 = zeros(1,128);
        pattern_order2 = zeros(1,128);
%     end
else
    pattern_order1 = zeros(1,120);
    pattern_order2 = zeros(1,120);
end
    
e_group_order = [1 5 2 6 3 7 4 8];
for ii = 1:8
    iGroup = e_group_order(ii);
    pattern_order1(ii:8:end) = e_groups{iGroup}(randperm(length(e_groups{iGroup})));
    pattern_order2(ii:8:end) = e_groups{iGroup}(randperm(length(e_groups{iGroup})));
end


%%
% first two chunks: single bipolar polarity combination

if elec_spacing == 60
    for ii = 1:32
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{1}(ii) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{1}(ii+32) = pInd;
        end
    end
    
    for ii = 33:64
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{2}(ii-32) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{2}(ii) = pInd;
        end
    end
    
    for ii = 65:96
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{3}(ii-64) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{3}(ii-32) = pInd;
        end
    end
    
    for ii = 97:128
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{4}(ii-96) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{4}(ii-64) = pInd;
        end
    end
    
    if both_polarity_combs

        for ii = 1:32
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{5}(ii) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{5}(ii+32) = pInd;
            end
        end
        
        for ii = 33:64
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{6}(ii-32) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{6}(ii) = pInd;
            end
        end
        
        for ii = 65:96
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{7}(ii-64) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{7}(ii-32) = pInd;
            end
        end
        for ii = 97:128
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{8}(ii-96) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{8}(ii-64) = pInd;
            end
        end
    end
    
    if include_single_elecs
        
        for ii = 1:32
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{5+4*both_polarity_combs}(ii) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{5+4*both_polarity_combs}(ii+32) = pInd;
            end
        end
        
        for ii = 33:64
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{6+4*both_polarity_combs}(ii-32) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{6+4*both_polarity_combs}(ii) = pInd;
            end
        end
        
        for ii = 65:96
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{7+4*both_polarity_combs}(ii-64) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{7+4*both_polarity_combs}(ii-32) = pInd;
            end
        end
        
        for ii = 97:128
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{8+4*both_polarity_combs}(ii-96) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{8+4*both_polarity_combs}(ii-64) = pInd;
            end
        end
    end
else
    for ii = 1:30
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{1}(ii) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{1}(ii+30) = pInd;
        end
    end
    
    for ii = 31:60
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{2}(ii-30) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{2}(ii) = pInd;
        end
    end
    
    for ii = 61:90
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{3}(ii-60) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{3}(ii-30) = pInd;
        end
    end
    
    for ii = 91:120
        pInd = find(elec_1==pattern_order1(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{4}(ii-90) = pInd;
        end
        pInd = find(elec_1==pattern_order2(ii) & stim_type==1);
        if ~isempty(pInd)
            pattern_order_all{4}(ii-60) = pInd;
        end
    end
    
    if both_polarity_combs
        for ii = 1:30
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{5}(ii) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{5}(ii+30) = pInd;
            end
        end
        
        for ii = 31:60
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{6}(ii-30) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{6}(ii) = pInd;
            end
        end
        
        for ii = 61:90
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{7}(ii-60) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{7}(ii-30) = pInd;
            end
        end
        for ii = 91:120
            pInd = find(elec_1==pattern_order1(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{8}(ii-90) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==2);
            if ~isempty(pInd)
                pattern_order_all{8}(ii-60) = pInd;
            end
        end
    end
    
    if include_single_elecs
        
        for ii = 1:30
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{5+4*both_polarity_combs}(ii) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{5+4*both_polarity_combs}(ii+30) = pInd;
            end
        end
        
        for ii = 31:60
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{6+4*both_polarity_combs}(ii-30) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{6+4*both_polarity_combs}(ii) = pInd;
            end
        end
        
        for ii = 61:90
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{7+4*both_polarity_combs}(ii-60) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{7+4*both_polarity_combs}(ii-30) = pInd;
            end
        end
        
        for ii = 91:120
            pInd = find(elec_1==pattern_order1(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{8+4*both_polarity_combs}(ii-90) = pInd;
            end
            pInd = find(elec_1==pattern_order2(ii) & stim_type==3);
            if ~isempty(pInd)
                pattern_order_all{8+4*both_polarity_combs}(ii-60) = pInd;
            end
        end
    end
end


%% 

% for jj = 1:length(pattern_order_all)
%     pattern_order_all{jj}(pattern_order_all{jj}==0) = [];
% end


%% check that each stimulus is applied twice
times_applied = zeros(1,size(array,2));
for ii = 1:length(pattern_order_all)
    for jj =1:length(pattern_order_all{ii})
        if pattern_order_all{ii}(jj)~=0
            times_applied(pattern_order_all{ii}(jj)) = times_applied(pattern_order_all{ii}(jj))+1;
        end
    end
end

if ~all(times_applied >= 2 ) %bad hack!
    error('each pattern should be applied twice...so something''g wrong!')
end

%%

%%

MovieChunks = [];
for ii = 1:length(pattern_order_all)
    thisPattern = pattern_order_all{ii};
    Times = (0:delay:delay*(length(thisPattern)-1)) + timeShift;
    
    %get rid of 0 patterns
    for jj = 1:length(thisPattern)
        Times(thisPattern == 0) = [];
        thisPattern(thisPattern == 0) = [];
    end
    
    if any(Times > nSamples)
        error('stimulus time exceeds movie chunk time')
    end
        
    Chunk=NS_MovieChunkGenerationForExperiment(Times,nSamples,thisPattern);
    MovieChunks = [MovieChunks Chunk];
end
MovieChunksFile=[length(pattern_order_all) MovieChunks]; %only one movie


if same_polarity
    figure; subplot(1,2,1); imagesc(array); title('opp polarity pairs')
    array(find(array == -1)) = 1; 
    subplot(1,2,2); imagesc(array); title('same polarity pairs');
end

if use_ratios
    figure; subplot(1,2,1); imagesc(array); title('opp polarity pairs')
    if scale_factor<0
        array(find(array == 1)) = -scale_factor;
    else
        array(find(array == -1)) = -scale_factor; 
    end
    subplot(1,2,2); imagesc(array); title('ratio');
end
arrayarray = cat(2,array1,array); 
figure; imagesc(arrayarray); 
% keyboard; 
%% Save files
% cd /Users/grosberg/Desktop/; 
allThresholds = ones(1,512); 
% allThresholds = load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/data_files/Safe1029 (1).mat');
% allThresholds = load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/data_files/Safe10296.mat');
% allThresholds = load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/data_files/Safe1193.mat');
allThresholds = load('/Users/grosberg/Downloads/Safe2016154.mat')
try
    allThresholds = allThresholds.Safe1029;
catch
    try
        allThresholds = allThresholds.Safe10296;
    catch
        try
            allThresholds = allThresholds.Safe1193;
        catch
            allThresholds_temp = ones(1,512); 
             allThresholds_temp(electrodes) = allThresholds.Safe2016154; 
             allThresholds = allThresholds_temp; 
        end
    end
end
meanThresh = mean(allThresholds(find(allThresholds))); 
 
allThresholds(allThresholds == 0) = meanThresh; 
if size(allThresholds,1) < size(allThresholds,2)
    allThresholds = allThresholds'; 
end
bundleThresholds = allThresholds(electrodes); 

%%

saveFiles = 0;
% saveName = 'safeZoneRatios'; 
list1 = [1    1   1    0.75 0.5 0.25 1]; 
list2 = [0.25 0.5 0.75 1    1   1    1]; 
figure; 
for j = 1:length(list1)
    arrayarrayhey =  arrayarray; 
    arrayarrayhey(find(arrayarrayhey == 1)) = list1(j);
    arrayarrayhey(find(arrayarrayhey == -1)) = list2(j);
    bundleThresholdsM = repmat(bundleThresholds,1,size(arrayarray,2));
    arrayScaled = arrayarrayhey.*bundleThresholdsM;
    %     figure; imagesc(arrayScaled); title(['ratio ' num2str(j)])
    hold on;
    scatter(arrayScaled(1,:),arrayScaled(5,:));
    if saveFiles
        fid = fopen([saveName '_ratio' num2str(j) '_pt'],'wb','ieee-le.l64');
        fwrite(fid,arrayScaled,'double');
        fclose(fid);
    end
end

%% Simulate the stimulus

figure; 
for pat = 1:size(arrayScaled,2); 
    cla;
    scatter(elec_coords(:,1),elec_coords(:,2),'k'); axis image; 
    e = find(arrayScaled(:,pat)); 
    hold on; 
    scatter(elec_coords(electrodes(e),1),elec_coords(electrodes(e),2),100,...
        arrayScaled(e,pat),'filled'); axis image;
    caxis([0 2.5]); 
    title(num2str(pat)); 
    pause; 
end