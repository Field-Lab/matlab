datarun=load_data('2012-08-09-3/data002');
datarun=load_params(datarun);

on_idx=datarun.vision.cell_types{1,1}.cell_ids;
off_idx=datarun.vision.cell_types{1,2}.cell_ids;

% Average radius of on cells
RF_area=0;
red=zeros(30,1);
blue=zeros(30,1);
green=zeros(30,1);
for i=on_idx
    cell_idx=find(datarun.cell_ids==i);
    rf=datarun.vision.sta_fits{cell_idx}.sd;
    RF_area=RF_area+rf(1)*rf(2);
    red=red+datarun.vision.timecourses(cell_idx).r;
    blue=blue+datarun.vision.timecourses(cell_idx).b;
    green=green+datarun.vision.timecourses(cell_idx).g; 
end
ON.RF_radius=sqrt(RF_area/length(on_idx));
ON.red=red/length(on_idx);
ON.blue=blue/length(on_idx);
ON.green=green/length(on_idx);

% Average radius of off cells
RF_area=0;
red=zeros(30,1);
blue=zeros(30,1);
green=zeros(30,1);
for i=off_idx
    cell_idx=find(datarun.cell_ids==i);
    rf=datarun.vision.sta_fits{cell_idx}.sd;
    RF_area=RF_area+rf(1)*rf(2);
    red=red+datarun.vision.timecourses(cell_idx).r;
    blue=blue+datarun.vision.timecourses(cell_idx).b;
    green=green+datarun.vision.timecourses(cell_idx).g; 
end
OFF.RF_radius=sqrt(RF_area/length(off_idx));
OFF.red=red/length(off_idx);
OFF.blue=blue/length(off_idx);
OFF.green=green/length(off_idx);

save('ON_STA_params.mat',ON);
save('OFF_STA_params.mat',OFF);
