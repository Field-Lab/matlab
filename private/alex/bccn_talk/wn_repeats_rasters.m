
%% 2015-08-17-1

% NDF4
% data002 BW-20-8-0.48-11111 1400s
% data003 BW-20-8-0.48-11111 repeats
% data004 BW-16-8-0.48-11111 repeats
% data005 NSEM repeats
% data006 BW-16-8-0.48-11111 1200s

% NDF3
% data007 BW-16-8-0.48-11111 1800s
% data008 BW-16-8-0.48-11111 repeats
% data009 BW-10-6-0.48-11111 repeats
% data010 NSEM repeats

% NDF2
% data011 BW-10-6-0.48-11111 1800s
% data012 BW-10-6-0.48-11111 repeats
% data013 BW-10-4-0.48-11111 repeats
% data014 NSEM repeats


%%
path2data = '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/';

% ndf3 and 2
for_use = ['data009'; 'data012'];

datarun = cell(2,1);
for i=1:2
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_neurons(datarun{i}); 
end

starun = load_data(fullfile(path2data, 'data002', 'data002'));
starun = load_params(starun,'verbose',1);


%%
cell_type = 4;
cellIDs = get_cell_indices(starun, {cell_type});

for k=1:length(cellIDs)
    
    figure
    set(gcf,'Name',[starun.cell_types{cell_type}.name, '  cell ID ', int2str(cellIDs(k))])
    set(gcf, 'position', [1   823   804   275]);
    for i=1:2
        spikes=datarun{i}.spikes{cellIDs(k)};
        trigs = datarun{i}.triggers;
        myTrigs=[0 find(diff(trigs)>0.84)'];
        splitSpikes=cell(14,1);
        for j=6:19
            tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
                - trigs(myTrigs(j)+1);
            tmp(tmp>10)=[];
            splitSpikes{j-5}=tmp*1000;
        end
        hold on
        for j=1:14
            plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13 + (i-1), '.k', 'markersize',0.1)
        end
    end
end

%%
path2data = '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/';

% ndf4 and 3
for_use = ['data004'; 'data009'];

datarun = cell(2,1);
for i=1:2
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_neurons(datarun{i}); 
end

starun = load_data(fullfile(path2data, 'data002', 'data002'));
starun = load_params(starun,'verbose',1);


%%
cell_type = 3;
cellIDs = get_cell_indices(starun, {cell_type});

for k=1:length(cellIDs)
    
    figure
    set(gcf,'Name',[starun.cell_types{cell_type}.name, '  cell ID ', int2str(cellIDs(k))])
    set(gcf, 'position', [1   823   804   275]);
    for i=1:2
        spikes=datarun{i}.spikes{cellIDs(k)};
        trigs = datarun{i}.triggers;
        myTrigs=[0 find(diff(trigs)>0.84)'];
        splitSpikes=cell(14,1);
        for j=6:19
            tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
                - trigs(myTrigs(j)+1);
            tmp(tmp>10)=[];
            splitSpikes{j-5}=tmp*1000;
        end
        hold on
        for j=1:14
            plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13 + (i-1), '.k', 'markersize',0.1)
        end
    end
end





%% 2015-03-09-2

% NDF2
% data015 BW-16-6-0.48-11111 900s
% data016 NSEM repeats
% data017 BW-8-6-0.48-11111 repeats
% data018 BW-8-6-0.48-11111 1200s

% NDF1
% data019 RGB-16-6-0.48-11111 900s
% data020 NSEM repeats
% data021 BW-8-6-0.48-11111 repeats
% data022 BW-8-6-0.48-11111 1200s


%%
path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/';

% ndf2 and 1
for_use = ['data017'; 'data021'];

datarun = cell(2,1);
for i=1:2
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_neurons(datarun{i}); 
end

starun = load_data(fullfile(path2data, 'data015', 'data015'));
starun = load_params(starun,'verbose',1);


