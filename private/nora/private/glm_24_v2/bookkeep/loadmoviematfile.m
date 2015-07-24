% AKHeitman 2014-03-03
% edit as we go along.. but nice to have a central loading file
% Clean short easy way to load the single mat file which usually
% contains the movies files in cell format (per block)
% and saved as uint8 ( may be different for BW)

% Output could be on a (0 1) scale or a (0 255) uint8 scale
%%%%%%%%%%

% updated 2014-04-06 .. added BD = NSEM_BaseDirectories, 
%                   .. added exp_info = experimentinfoNSEM(exp_nm)

% test code
%{
exp_nm = '2012-08-09-3'; stim_type = 'NSEM'; 
cone_model = '8pix_Model1_1e4_8pix'; movietype = 'testmovie';
[blockedmoviecell, fileloaded] = loadmoviematfile(exp_nm , stim_type, cone_model,movietype);
%}
function [blockedmoviecell,  novelmoviestats, origmatfile] = loadmoviematfile(exp_nm , stim_type, cone_model,movietype)

BD = NSEM_BaseDirectories;
stimuli_basedir = BD.NSEM_stimuli;
[exp_info] = experimentinfoNSEM(exp_nm);

if strcmp(stim_type, 'WN') || strcmp(stim_type, 'BW')
    moviename = exp_info.WNfile; scheme = '';
end
if strcmp(stim_type, 'NSEM')
    moviename = exp_info.NSEMfile; scheme = exp_info.NSEMscheme;
end

%{
if strcmp(stim_type,'BW') || strcmp(stim_type, 'WN')
    moviename = 'BW-8-1-0.48-11111_RNG_16807';
    scheme = '';
end

if strcmp(stim_type, 'NSEM')    
    if strcmp(exp_nm, '2012-08-09-3') ||  strcmp(exp_nm,'2012-09-27-3')
        moviename = 'NSEM_eye-120-3_0-3600';
        scheme = 'schemeA_';
    end
    
    if strcmp(exp_nm, '2013-08-19-6')
        moviename = 'NSEM_eye-long-v2';
        scheme = 'schemeA_';
    end
    if strcmp(exp_nm, '2013-10-10-0')
        moviename = 'NSEM_FEM900FF_longrast';
        scheme = 'schemeB_';
    end
end
%}


% fit : BWmovie  NSEMmovie  ( BWmovie.fitmovie.movie_byblock)
% test: testmovie.matrix

if strcmp(stim_type, 'WN')
    origmatfile = sprintf('%s/%s/%s_%s%s.mat' , stimuli_basedir, moviename, movietype , scheme , cone_model);
end
if strcmp(stim_type, 'NSEM')
    origmatfile = sprintf('%s/%s/%s_%s_%s.mat' , stimuli_basedir, moviename, movietype , scheme , cone_model);
end
    


eval(sprintf('load %s' ,origmatfile ))

if strcmp(movietype, 'fitmovie')
    if strcmp(stim_type, 'NSEM')
        blockedmoviecell = NSEMmovie.fitmovie.movie_byblock;
    elseif strcmp(stim_type, 'BW') || strcmp(stim_type, 'WN')
        blockedmoviecell = BWmovie.fitmovie.movie_byblock;
    end
end
if strcmp(movietype, 'testmovie')
    blockedmoviecell{1} = testmovie;
end


if strcmp(stim_type, 'NSEM')
    statsfile = sprintf('%s/%s/inputstats_%s.mat' , stimuli_basedir, moviename, cone_model);
    eval(sprintf('load %s' ,statsfile ));
    novelmoviestats = inputstats;
end
if strcmp(stim_type, 'BW') || strcmp(stim_type, 'WN')
    novelmoviestats.mu_avgIperpix = .5;
end


%{
if strcmp(stim_type, 'BW') && strcmp(cone_model, '8pix_Identity_8pix')
        load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/BW-8-1-0.48-11111_RNG_16807/BWmovie.mat
        fullmovie_struct = BWmovie; clear BWmovie
end

if strcmp(stim_type, 'BW') && strcmp(cone_model, '8pix_Model1_1e4_8pix')
        load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/BW-8-1-0.48-11111_RNG_16807/fitmovie_8pix_Model1_1e4_8pix.mat
        fullmovie_struct = BWmovie; clear BWmovie
        fullmovie_struct.params.fit_frames = 3600;
        fullmovie_struct.fitmovie.ind_to_block = 2:2:120;
end

if strcmp(cone_model, '8pix_Identity_8pix') || strcmp(cone_model, '8pix_PtLog_8pix')
    if strcmp(stim_type, 'NSEM')
        if strcmp(exp_nm, '2012-08-09-3') || strcmp(exp_nm,'2012-09-27-3') || strcmp(exp_nm, '2012-04-13-4')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600/inputstats_8pix_Identity_8pix.mat
        end

        if strcmp(exp_nm, '2013-08-19-6')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-long-v2/fitmovie_schemeA_8pix_Identity_8pix.mat
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-long-v2/inputstats_8pix_Identity_8pix.mat
        end


        if strcmp(exp_nm, '2013-10-10-0')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_FEM900FF_longrast/fitmovie_schemeB_8pix_Identity_8pix.mat
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_FEM900FF_longrast/inputstats_8pix_Identity_8pix.mat
        end

        fullmovie_struct = NSEMmovie;    clear NSEMmovie;
    end
    if strcmp( exp_nm , '2013-10-10-0')
        fullmovie_struct.params.width = 80 ;
        fullmovie_struct.params.height = 40;
        
        if strcmp(fit_type, 'NSEM')
            fullmovie_struct.params.fit_frames = 14400;
        end
        if strcmp(fit_type, 'BW')
            fullmovie_struct.params.fit_frames = 3600;
        end
            
            
    end
end

if strcmp(cone_model, '8pix_Model1_1e4_8pix') && strcmp(stim_type ,  'NSEM')
    if strcmp(exp_nm, '2012-08-09-3') || strcmp(exp_nm,'2012-09-27-3') || strcmp(exp_nm, '2012-04-13-4')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Model1_1e4_8pix.mat 
             load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600/inputstats_8pix_Model1_1e4_8pix.mat
            fullmovie_struct = NSEMmovie; 
            fullmovie_struct.params.width = 80 ;
            fullmovie_struct.params.height = 40;
            fullmovie_struct.params.fit_frames = 7200;
    end
    if strcmp(exp_nm, '2013-08-19-6')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-long-v2/fitmovie_schemeA_8pix_Model1_1e4_8pix.mat 
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-long-v2/inputstats_8pix_Model1_1e4_8pix.mat
           
            fullmovie_struct = NSEMmovie; 
            fullmovie_struct.params.width = 80 ;
            fullmovie_struct.params.height = 40;
            fullmovie_struct.params.fit_frames = 7200;
    end
    
    if strcmp(exp_nm, '2013-10-10-0')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_FEM900FF_longrast/fitmovie_schemeB_8pix_Model1_1e4_8pix.mat 
             load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_FEM900FF_longrast/inputstats_8pix_Model1_1e4_8pix.mat
            fullmovie_struct = NSEMmovie; 
            fullmovie_struct.params.width = 80 ;
            fullmovie_struct.params.height = 40;
            fullmovie_struct.params.fit_frames = 14400;
    end
    
end

if strcmp(cone_model, '8pix_Model1_1e5_8pix') && strcmp(stim_type ,  'NSEM')
    if strcmp(exp_nm, '2012-08-09-3') || strcmp(exp_nm,'2012-09-27-3') || strcmp(exp_nm, '2012-04-13-4')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Model1_1e5_8pix.mat 
             load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600/inputstats_8pix_Model1_1e5_8pix.mat
            fullmovie_struct = NSEMmovie; 
            fullmovie_struct.params.width = 80 ;
            fullmovie_struct.params.height = 40;
            fullmovie_struct.params.fit_frames = 7200;
    end
    if strcmp(exp_nm, '2013-08-19-6')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-long-v2/fitmovie_schemeA_8pix_Model1_1e5_8pix.mat 
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-long-v2/inputstats_8pix_Model1_1e5_8pix.mat
           
            fullmovie_struct = NSEMmovie; 
            fullmovie_struct.params.width = 80 ;
            fullmovie_struct.params.height = 40;
            fullmovie_struct.params.fit_frames = 7200;
    end
    
    if strcmp(exp_nm, '2013-10-10-0')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_FEM900FF_longrast/fitmovie_schemeB_8pix_Model1_1e5_8pix.mat 
             load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_FEM900FF_longrast/inputstats_8pix_Model1_1e5_8pix.mat
            fullmovie_struct = NSEMmovie; 
            fullmovie_struct.params.width = 80 ;
            fullmovie_struct.params.height = 40;
            fullmovie_struct.params.fit_frames = 14400;
    end
    
end


if strcmp(cone_model, '8pix_Model1_1e6_8pix') && strcmp(stim_type ,  'NSEM')
    if strcmp(exp_nm, '2012-08-09-3') || strcmp(exp_nm,'2012-09-27-3') || strcmp(exp_nm, '2012-04-13-4')
            load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Model1_1e6_8pix.mat 
            fullmovie_struct = NSEMmovie; 
            fullmovie_struct.params.width = 80 ;
            fullmovie_struct.params.height = 40;
            fullmovie_struct.params.fit_frames = 7200;
    end   
end
%}
end