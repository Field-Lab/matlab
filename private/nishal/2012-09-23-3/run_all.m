%%
load('pattern_cell_dat.csv');


patt_list = pattern_cell_dat(:,1);
cell_list = pattern_cell_dat(:,2);

for istim = [5:length(patt_list),1]
    fname=sprintf('/Volumes/Analysis/nishal/2012-09-24-3/stim%d',istim);
    
    if(exist(fname)~=7)
    mkdir(fname);
    end
    
patternNo = patt_list(istim);
cell_no = cell_list(istim);

find_data2;
end