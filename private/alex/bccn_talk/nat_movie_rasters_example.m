clear

%% 2015-03-09-2
% NDF5
% data005 BW-20-10-0.48-11111-16x16 1800s
% data006 NSEM repeats
% data007 BW-20-10-0.48-11111 repeats

% NDF4
% data008 BW-10-8-0.48-11111 1800s
% data009 NSEM repeats
% data010 BW-10-8-0.48-11111 repeats

% NDF3
% data011 BW-16-8-0.48-11111 1800s
% data012 NSEM repeats
% data013 BW-8-8-0.48-11111 repeats
% data014 BW-8-8-0.48-11111 1386s

% NDF2
% data015 BW-16?6-0.48-11111 900s
% data016 NSEM repeats
% data017 BW-8-6-0.48-11111 repeats
% data018 BW-8-6-0.48-11111 1200s

% NDF1
% data019 RGB-16-6-0.48-11111 900s
% data020 NSEM repeats
% data021 BW-8-6-0.48-11111 repeats
% data022 BW-8-6-0.48-11111 1200s

% NDF0
% data023 RGB-16-4-0.48-11111 900s
% data024 NSEM repeats
% data025 BW-4-2-0.48-11111 repeats
% data026 BW-4-2-0.48-11111 1800s
% data027 RGB-4-2-0.48-11111 1800s

%%
% load('/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/data_011_015_chopped.mat','data1','data2')

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/';


for_use = ['data006'; 'data009'; 'data012'; 'data016'; 'data020'; 'data024'];

datarun = cell(6,1);
for i=1:6
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_neurons(datarun{i}); 
end

starun = load_data(fullfile(path2data, 'data005', 'data005'));
starun = load_params(starun,'verbose',1);

%
% nds={'NDF3', 'NDF2'};
% sr=120;
% sig=50;
% st=10000/sr*6.2/(60*sig);
% time_list=-3.1:st:3.1;
% kern=zeros(1,length(time_list));
% for i=1:length(time_list)
%     kern(i)=250/sig*exp((1-time_list(i)^2)/2);
% end
%
% figure
% plot(kern)



%% OFF midget, cell 2945

% cellID = 2945;
% datarunID=find(starun.cell_ids==cellID);

cell_type = 1;
cellIDs = get_cell_indices(starun, {cell_type});

for k=1:length(cellIDs)
    
    figure
    set(gcf,'Name',[starun.cell_types{cell_type}.name, '  cell ID ', int2str(cellIDs(k))])
    
    for i=2:6
        spikes=datarun{i}.spikes{cellIDs(k)};
        trigs = datarun{i}.triggers;
        myTrigs=[0 find(diff(trigs)>0.9)'];
        splitSpikes=cell(14,1);
        for j=6:19
            tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
                - trigs(myTrigs(j)+1);
            splitSpikes{j-5}=tmp*1000;
        end
        hold on
        for j=1:14
            plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13 + (i-1), '.k', 'markersize',0.1)
        end
    end
end

onp = [1 2 5 8 9 11 13 16 18 22 24 25 27 28];
offp = [2 4 6 9 10 11 12 14 16 18 19 24 25 28 29 30 31 32];
onm = [88 83 79 76 75 73 65 63 54 53 51 50 48 45 43 34 30 29];
offm = [2 4 6 8 9 16 28 30 32 34 38 41];




