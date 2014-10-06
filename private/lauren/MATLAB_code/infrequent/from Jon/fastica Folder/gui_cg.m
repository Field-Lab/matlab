function gui_cg()

% This function is needed by FASTICAG

% This file just removes the global variables
% that are used in FASTICAG from the memory

% 2.4.1998
% Hugo Gävert

clear global c_FastICA_appr_strD;
clear global c_FastICA_appr_strV;
clear global c_FastICA_dMod_strD;
clear global c_FastICA_dMod_strV;
clear global c_FastICA_g1_strD;
clear global c_FastICA_g1_strV;
clear global c_FastICA_g2_strD;
clear global c_FastICA_g2_strV;
clear global c_FastICA_iSta_strD;
clear global c_FastICA_iSta_strV;
clear global c_FastICA_verb_strD;
clear global c_FastICA_verb_strV;
clear global g_FastICA_a1;
clear global g_FastICA_a2;
clear global g_FastICA_approach;
clear global g_FastICA_displayIn;
clear global g_FastICA_displayMo;
clear global g_FastICA_epsilon;
clear global g_FastICA_g;
clear global g_FastICA_ica_A;
clear global g_FastICA_ica_W;
clear global g_FastICA_ica_sig;
clear global g_FastICA_initGuess;
clear global g_FastICA_initState;
clear global g_FastICA_loadType;
clear global g_FastICA_maxNumIte;
clear global g_FastICA_mixedsig;
clear global g_FastICA_numOfIC;
clear global g_FastICA_pca_D;
clear global g_FastICA_pca_E;
clear global g_FastICA_verbose;
clear global g_FastICA_white_dwm;
clear global g_FastICA_white_sig;
clear global g_FastICA_white_wm;
clear global g_mix_meanremoved;
clear global h_b_advOpt;
clear global h_b_initGuess;
clear global h_b_loaddata;
clear global h_e_a1;
clear global h_e_a2;
clear global h_e_displayInterval;
clear global h_e_epsilon;
clear global h_e_file;
clear global h_e_suffix;
clear global h_e_maxIterations;
clear global h_e_numOfIC;
clear global h_f_FastICA;
clear global h_f_FastICAAdvOpt;
clear global h_f_FastICALoad;
clear global h_f_FastICASave;
clear global h_ica_A;
clear global h_ica_W;
clear global h_ica_sig;
clear global h_initGuess;
clear global h_mixed;
clear global h_pca_D;
clear global h_pca_E;
clear global h_pm_approach;
clear global h_pm_displayMode;
clear global h_pm_g;
clear global h_pm_initState;
clear global h_pm_verbose;
clear global h_t_dim;
clear global h_t_icaStatus;
clear global h_t_initGuess;
clear global h_t_mixedStatus;
clear global h_t_newDim;
clear global h_t_newdim;
clear global h_t_numOfSamp;
clear global h_t_numofsamp;
clear global h_t_whiteStatus;
clear global h_white_dwm;
clear global h_white_sig;
clear global h_white_wm;