% load fitted results

path = '/Volumes/Lab/Users/crhoades/Fitting/Test/2010-09-24-0/data001-nwpca/ON midget';

files=dir(fullfile(path,'*.mat')); 

for i = 1:size(files,1)
%     if isempty(strfind(files(i).name, '7687')) && isempty(strfind(files(i).name, '5463'))
    load([path, '/', files(i).name])
    
    rf_diameter_center{1}(i) = max(parameters_per_cell.center_sd_x, parameters_per_cell.center_sd_y)*640/parameters_per_cell.x_dim*2*5;
    surround_strength{1}(i) = parameters_per_cell.surround_amp_scale;
%     end
    
end


path = '/Volumes/Lab/Users/crhoades/Fitting/2010-09-24-0/data001-nwpca/OFF midget';

files=dir(fullfile(path,'*.mat')); 

for i = 1:size(files,1)
    if isempty(strfind(files(i).name, '1636')) && isempty(strfind(files(i).name, '2236')) && isempty(strfind(files(i).name, '2312')) && isempty(strfind(files(i).name, '2762')) && isempty(strfind(files(i).name, '6482')) && isempty(strfind(files(i).name, '7112')) && isempty(strfind(files(i).name, '4516')) && isempty(strfind(files(i).name, '5207')) 
    load([path, '/', files(i).name])
    
    rf_diameter_center{2}(i) = max(parameters_per_cell.center_sd_x, parameters_per_cell.center_sd_y)*640/parameters_per_cell.x_dim*2*5;
    surround_strength{2}(i) = parameters_per_cell.surround_amp_scale;
    end
    
end

figure; 
[nb,xb] = hist(rf_diameter_center{1}); 
b = bar(xb,nb);
hold on; 
[nb,xb] =hist(rf_diameter_center{2});
bh= bar(xb,nb);
set(bh, 'facecolor', [0 1 0])
legend('ON midget', 'OFF midget')
