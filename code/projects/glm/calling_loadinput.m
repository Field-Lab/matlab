clear
exps = 3;
stimtypes = [1]; % white noise only  (2 is natural scens)
celltypes = [1]; % only ON Parasol
cell_subset = 'debug';
glm_settings{1}.type = 'debug';
glm_settings{1}.name = 'true';
glm_settings{2}.type = 'cone_model';
glm_settings{2}.name = 'rieke_linear';
runoptions.load_previousfit.glm_settings = glm_settings;
glm_settings{3}.type= 'input_pt_nonlinearity';
glm_settings{3}.name= 'log_powerraise';
glm_wrap(exps,stimtypes,celltypes,cell_subset,glm_settings,{},runoptions)



%%%
clear
exps = 3;
stimtypes = [1]; % white noise only  (2 is natural scens)
celltypes = [1]; % only ON Parasol
cell_subset = 'debug';
glm_settings{1}.type = 'debug';
glm_settings{1}.name = 'true';
runoptions.replace_existing = true;
glm_wrap(exps,stimtypes,celltypes,cell_subset,glm_settings,{},runoptions)
%%% Should have the following minimization sequence  
### running: WN expC ONPar_2824: debug_fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams ###

                                Norm of      First-order 
 Iteration        f(x)          step          optimality   CG-iterations
     0            1297.56                         2e+04                
     1            1297.56             10          2e+04           4
     2           -42251.5            2.5       2.56e+03           0
     3           -45320.3        4.13912       5.15e+03           7