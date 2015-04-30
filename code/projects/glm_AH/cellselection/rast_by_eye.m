% AKHeitman 2014-12-10
% Hack Code to Transfer Good/Bad cells example from SPreadsheet
% Into a format that is readily accesible


clear; clc

BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
rast_eye = allcells;


eval(sprintf('load %s/Raster_By_Eye/rastcheck_byeye_expA.mat',BD.Cell_Selection));
eval(sprintf('load %s/Raster_By_Eye/rastcheck_byeye_expB.mat',BD.Cell_Selection));
eval(sprintf('load %s/Raster_By_Eye/rastcheck_byeye_expC.mat',BD.Cell_Selection));
eval(sprintf('load %s/Raster_By_Eye/rastcheck_byeye_expD.mat',BD.Cell_Selection));


%%

 
rast_eye{1}.stim_type{1}.celltype{1}.not_stable = A_WN_ONP(:,2); 
rast_eye{1}.stim_type{1}.celltype{2}.not_stable = A_WN_OFFP(:,2); 
rast_eye{1}.stim_type{2}.celltype{1}.not_stable = A_NSEM_ONP(:,2); 
rast_eye{1}.stim_type{2}.celltype{2}.not_stable = A_NSEM_OFFP(:,2); 
rast_eye{1}.stim_type{1}.celltype{1}.not_sharp = A_WN_ONP(:,3); 
rast_eye{1}.stim_type{1}.celltype{2}.not_sharp = A_WN_OFFP(:,3); 
rast_eye{1}.stim_type{2}.celltype{1}.not_sharp = A_NSEM_ONP(:,3); 
rast_eye{1}.stim_type{2}.celltype{2}.not_sharp = A_NSEM_OFFP(:,3); 
rast_eye{1}.stim_type{1}.celltype{1}.sharp = A_WN_ONP(:,4); 
rast_eye{1}.stim_type{1}.celltype{2}.sharp = A_WN_OFFP(:,4); 
rast_eye{1}.stim_type{2}.celltype{1}.sharp = A_NSEM_ONP(:,4); 
rast_eye{1}.stim_type{2}.celltype{2}.sharp = A_NSEM_OFFP(:,4); 

rast_eye{2}.stim_type{1}.celltype{1}.not_stable = B_WN_ONP(:,2); 
rast_eye{2}.stim_type{1}.celltype{2}.not_stable = B_WN_OFFP(:,2); 
rast_eye{2}.stim_type{2}.celltype{1}.not_stable = B_NSEM_ONP(:,2); 
rast_eye{2}.stim_type{2}.celltype{2}.not_stable = B_NSEM_OFFP(:,2); 
rast_eye{2}.stim_type{1}.celltype{1}.not_sharp = B_WN_ONP(:,3); 
rast_eye{2}.stim_type{1}.celltype{2}.not_sharp = B_WN_OFFP(:,3); 
rast_eye{2}.stim_type{2}.celltype{1}.not_sharp = B_NSEM_ONP(:,3); 
rast_eye{2}.stim_type{2}.celltype{2}.not_sharp = B_NSEM_OFFP(:,3); 
rast_eye{2}.stim_type{1}.celltype{1}.sharp = B_WN_ONP(:,4); 
rast_eye{2}.stim_type{1}.celltype{2}.sharp = B_WN_OFFP(:,4); 
rast_eye{2}.stim_type{2}.celltype{1}.sharp = B_NSEM_ONP(:,4); 
rast_eye{2}.stim_type{2}.celltype{2}.sharp = B_NSEM_OFFP(:,4); 

rast_eye{3}.stim_type{1}.celltype{1}.not_stable = C_WN_ONP(:,2); 
rast_eye{3}.stim_type{1}.celltype{2}.not_stable = C_WN_OFFP(:,2); 
rast_eye{3}.stim_type{2}.celltype{1}.not_stable = C_NSEM_ONP(:,2); 
rast_eye{3}.stim_type{2}.celltype{2}.not_stable = C_NSEM_OFFP(:,2); 
rast_eye{3}.stim_type{1}.celltype{1}.not_sharp = C_WN_ONP(:,3); 
rast_eye{3}.stim_type{1}.celltype{2}.not_sharp = C_WN_OFFP(:,3); 
rast_eye{3}.stim_type{2}.celltype{1}.not_sharp = C_NSEM_ONP(:,3); 
rast_eye{3}.stim_type{2}.celltype{2}.not_sharp = C_NSEM_OFFP(:,3); 
rast_eye{3}.stim_type{1}.celltype{1}.sharp = C_WN_ONP(:,4); 
rast_eye{3}.stim_type{1}.celltype{2}.sharp = C_WN_OFFP(:,4); 
rast_eye{3}.stim_type{2}.celltype{1}.sharp = C_NSEM_ONP(:,4); 
rast_eye{3}.stim_type{2}.celltype{2}.sharp = C_NSEM_OFFP(:,4); 

rast_eye{4}.stim_type{1}.celltype{1}.not_stable = D_WN_ONP(:,2); 
rast_eye{4}.stim_type{1}.celltype{2}.not_stable = D_WN_OFFP(:,2); 
rast_eye{4}.stim_type{2}.celltype{1}.not_stable = D_NSEM_ONP(:,2); 
rast_eye{4}.stim_type{2}.celltype{2}.not_stable = D_NSEM_OFFP(:,2); 
rast_eye{4}.stim_type{1}.celltype{1}.not_sharp = D_WN_ONP(:,3); 
rast_eye{4}.stim_type{1}.celltype{2}.not_sharp = D_WN_OFFP(:,3); 
rast_eye{4}.stim_type{2}.celltype{1}.not_sharp = D_NSEM_ONP(:,3); 
rast_eye{4}.stim_type{2}.celltype{2}.not_sharp = D_NSEM_OFFP(:,3); 
rast_eye{4}.stim_type{1}.celltype{1}.sharp = D_WN_ONP(:,4); 
rast_eye{4}.stim_type{1}.celltype{2}.sharp = D_WN_OFFP(:,4); 
rast_eye{4}.stim_type{2}.celltype{1}.sharp = D_NSEM_ONP(:,4); 
rast_eye{4}.stim_type{2}.celltype{2}.sharp = D_NSEM_OFFP(:,4); 
%%

eval(sprintf('save %s/Raster_By_Eye/rast_eye.mat rast_eye',BD.Cell_Selection));




            
            
            