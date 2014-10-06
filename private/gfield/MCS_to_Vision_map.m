%% load a RRS map

datarun = load_data('/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data008/data008');
datarun = load_neurons(datarun);
datarun = load_ei(datarun, 4, 'array_type', 519);

figure(1); clf; hold on;
plot(datarun.ei.position(:,1), datarun.ei.position(:,2), 'go')
for chnl = 1:512
    text(datarun.ei.position(chnl,1), datarun.ei.position(chnl,2), num2str(chnl))
end



% make the RRS Electrode array square to conform to the MCS array
Y_vals = -450:120:450
temp_rrs_positions = datarun.ei.position;
for vl = 1:length(Y_vals)
    temp_indices = find(temp_rrs_positions(:,2) == Y_vals(vl));
    temp_rrs_positions(temp_indices,1) = temp_rrs_positions(temp_indices,1) + 30;
end

% plot new electrode positions to check
figure(3); clf; hold on;
plot(temp_rrs_positions(:,1), temp_rrs_positions(:,2), 'go')
for chnl = 1:512
    text(temp_rrs_positions(chnl,1), temp_rrs_positions(chnl,2), num2str(chnl))
end


% make the minimum y value 0 to conform to MCS array
y_shift = -1*min(temp_rrs_positions(:,2)) * ones(512,1);
x_shift = -1*min(temp_rrs_positions(:,1)) * ones(512,1);
shifted_rrs_positions = temp_rrs_positions;
shifted_rrs_positions(:,2) = shifted_rrs_positions(:,2) +  y_shift;
shifted_rrs_positions(:,1) = shifted_rrs_positions(:,1) +  x_shift;


% plot new electrode positions to check
figure(3); clf; hold on;
plot(shifted_rrs_positions(:,1), shifted_rrs_positions(:,2), 'go')
for chnl = 1:512
    text(shifted_rrs_positions(chnl,1), shifted_rrs_positions(chnl,2), num2str(chnl+1))
end

%% MCS map from Neuroshare API


% data reading path
%cd /Volumes/DiskU/2011-08-01-4/data001/
%cd /Data/MCS/gfield/2011-08-04-1/data002/
cd /Data/MCS/gfield/2011-08-18-0/data000/

% data writing path
%path_to_write = '/Volumes/Palace/Data/gfield/2011-08-01-4/data001/';
%file_to_write = 'data000000.bin';
filename{1} = 'data000000.mcd';

% set library path
[pathname, name, ext]=fileparts(which('nsMCDLibrary.dylib'));
ns_SetLibrary([pathname filesep name ext]);


% Open data file
[nsresult, hfile] = ns_OpenFile(filename{1});
[nsresult, FileInfo] = ns_GetFileInfo(hfile);


MCS_positions = zeros(252,2);
for chnl = 1:252;
    [nsresult, temp_el] = ns_GetAnalogInfo(hfile, chnl);

    MCS_positions(chnl, :) = [temp_el.LocationX, temp_el.LocationY];
end




conversion_factor = 960/0.0016;
MCS_positions = MCS_positions * conversion_factor;

figure; clf; hold on; 
plot(MCS_positions(:,1), MCS_positions(:,2), 'ro')
for chnl = 1:252
    text(MCS_positions(chnl,1), MCS_positions(chnl,2), num2str(chnl))
end
axis([-40 1000 -40 1000])


new_map(1:252,1) = 1:252;
new_map(:,2) = MCS_positions(:,1);
new_map(:,3) = MCS_positions(:,2);

[nsresult, ContinuousCount, data] = ns_GetAnalogData(hfile, 3, 1, 20000*10);


%% load the MCS map

cd ~/Desktop
load mcs_array_map


[~,sorted_indices] = sort(mcs_array_map(:,1), 'ascend');
mcs_array_map = mcs_array_map(sorted_indices,:);

figure(10); clf; hold on;
plot(mcs_array_map(:,2), mcs_array_map(:,3), 'ro')
for chnl = 1:252
    text(mcs_array_map(chnl,2), mcs_array_map(chnl,3), num2str(mcs_array_map(chnl,1)))
end


new_map = mcs_array_map;
new_map(:,3) = new_map(:,3) - 900;
new_map(:,3) = new_map(:,3) *-1;

figure(11); clf; hold on;
plot(new_map(:,2), new_map(:,3), 'ro')
for chnl = 1:252
    text(new_map(chnl,2), new_map(chnl,3), num2str(chnl+1))
end


%% Converted MCS support email

% set path
%[nsresult] = ns_SetLibrary('nsMCDlibrary.dll')
[pathname, name, ext]=fileparts(which('nsMCDLibrary.dylib'));
ns_SetLibrary([pathname filesep name ext]);

% open data file
cd /Data/MCS/gfield/2011-08-18-0/data000/
filename{1} = 'data000000.mcd';

[nsresult, hfile] = ns_OpenFile(filename{1})
[nsresult,info]=ns_GetFileInfo(hfile)

figure
for i = 1:info.EntityCount
    [nsresult,entity] = ns_GetEntityInfo(hfile,i);
    if ( entity.EntityType == 2 && strcmp(entity.EntityLabel(1:8),'elec0001'))
        [nsresult,analog]=ns_GetAnalogInfo(hfile,i);
        text(analog.LocationX * 100, analog.LocationY * 100, entity.EntityLabel(26:27));
    end
end


%%

% get distances between RRS and MCS electrodes
dist_mat = ipdm(shifted_rrs_positions, new_map(:,2:3));

array_correspondences = zeros(513,1);
noise_pointer = 1;
noise_index = 254:1:513;
% Draw correspondences between RRS and MCS maps
for chnl = 1:253
    
    if chnl == 1
        array_correspondences(chnl) = 1;
    else
    
        temp_ind = find(dist_mat(:,chnl-1) < 0.001);
    
        if isempty(temp_ind)
            array_correspondences(chnl) = noise_index(noise_pointer);
            noise_pointer = noise_pointer +1;
        else
            array_correspondences(chnl) = temp_ind+1; 
        end
    end
end

electrode_counter = 1:513;
unassigned_electrodes = setdiff(electrode_counter, array_correspondences);

array_correspondences(254:513) = unassigned_electrodes;
mapping = array_correspondences;

cd ~/MCS-conversion/
save mapping mapping


%% sanity check
figure(13); clf; hold on
for chnl = 2:513
    plot(shifted_rrs_positions(array_correspondences(chnl)-1,1), shifted_rrs_positions(array_correspondences(chnl)-1,2), 'ro')

    text(shifted_rrs_positions(array_correspondences(chnl)-1,1), shifted_rrs_positions(array_correspondences(chnl)-1,2), num2str(chnl))
end




%%



% get distances between RRS and MCS electrodes
dist_mat = ipdm(shifted_rrs_positions, new_map(:,2:3));

clear mapping
for chnl = 1:252
    mapping(chnl) = find(dist_mat(:,chnl) == 0);
end

mapping = mapping +1;

electrode_list = 2:513;
missing_electrodes = setdiff(electrode_list, mapping);

% append extra electrodes that don't have a real correspondence
mapping = [mapping, missing_electrodes];

% prepend TTL channel
mapping = [1, mapping];











