% Version Convergence -- AKHeitman 2014-10-19

% Version 22-- AKHeitman 2014-04-02
%             -- done 2014-04-06
% Calls:  BD = NSEM_BaseDirectories
        % NSEM_secondaryDirectories
        %StimulusPars.master = params_classificationfile(exp_info.masfile);
        %StimulusPars.slv    = params_fitandraster(slv_file);


% MFILE     Directories_Params_func  ( so that we input date and do less
%                                                      switching)
%
%
% inputs: string_date... piece number in the char form '2012-09-27-3'
%         string_fitttyepe  .. either char 'BW'   or char 'NSEM'
%         boolean_debug   either false or true,  false = long version
%                             ture is the short debugger version
%         maptype:   either  "mapEI"  or "mapPRJ"
%
%

% slv_type is an optional argument
%
% usage:  single m-file with all directories and params for the GLM
%
%
% outputs:     datarun  -loads the master classifcation (1)
%                       -loads the slave classification and spikes (2)
%              all params used throughout the GLM  
%
% AK HEITMAN START 2012-09-17 .. COMMENTING 2012-10-1  


%%%%%%%%   NOTES      %%%%%%%%%%%
% DEBUG JUST LIMITS THE NUMBER OF BLOCKS
% NSEM AND BW HAVE DIFFERENT NUMBERING SCHEME IN MATLAB STORAGE
% NOVEL - ALL NOVEL BLOCKS
% TEST  - THE STATIC RUNS
% FIT   - SUBSET OF NOVEL USED FOR MODEL FITTING
% BW AND NSEM PARAMS BASED UPON MARTIN'S STIMULUS FILES

%{
clear all
string_piece = '2012-08-09-3'; slv_type = 'NSEM'; map_type = 'mapPRJ'; 
Directories_Params_v22(string_piece, slv_type, map_type)
[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v22(string_piece, slv_type, map_type)
%}
function [StimulusPars] = Directories_Params_v23_Conv(string_piece, map_type, percentage)

%%
%%%%% 0.   BASE DIRECTORY THAT NEEDS TO BE SET %%%%%%

[BD] = NSEM_BaseDirectories;
if ~exist('map_type', 'var'),  map_type = 'mapPRJ';  end
[exp_info] = experimentinfoNSEM(string_piece) ;
if strcmp(slv_type , 'WN') || strcmp(slv_type, 'BW')
    dn_slv = exp_info.dr.slvWN;
elseif strcmp(slv_type, 'NSEM')
    dn_slv = exp_info.dr.slvNSEM;
end
if strcmp(map_type, 'mapPRJ')
    dn_slv  = sprintf('%s-from-%s', dn_slv , exp_info.dr.mas);
end
clear inputs analysisdir directory_type 


%% find stimulus parameters
if strcmp(slv_type , 'WN') || strcmp(slv_type, 'BW')
    slv_file = exp_info.WNfile;
elseif strcmp(slv_type, 'NSEM')
    slv_file = exp_info.NSEMfile;
end
StimulusPars.master = stimparams_classificationfile(exp_info.masfile);
StimulusPars.slv    = stimparams_fitandraster(slv_file);
clear slv_file
%% experiment dependent fit/test parameters

StimulusPars.slv.FitBlocks       = 6:2:(StimulusPars.slv.n_blk-2);
StimulusPars.slv.TestBlocks      = 5:2:(StimulusPars.slv.n_blk-2);
if strcmp(string_piece, '2013-10-10-0') && strcmp(slv_type, 'NSEM')
    StimulusPars.slv.n_rep = 27;
    StimulusPars.slv.n_blk = 54;
    StimulusPars.slv.NovelBlocks     =  2:2:StimulusPars.slv.n_blk; 
    StimulusPars.slv.StaticBlocks    =  1:2:(StimulusPars.slv.n_blk-1);
    StimulusPars.slv.FitBlocks       = 4:2:54;
    StimulusPars.slv.TestBlocks      = 3:2:53;
end


old_length = length(StimulusPars.slv.FitBlocks);
new_length = max(floor( percentage*old_length ), 1);

StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:new_length); 




%%

end

