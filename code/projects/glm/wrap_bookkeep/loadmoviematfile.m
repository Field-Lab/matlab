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
    statsfile = sprintf('%s/%s/inputstats_%s.mat' , stimuli_basedir, moviename, cone_model);
    eval(sprintf('load %s' ,statsfile ));
    novelmoviestats = inputstats;
end


end