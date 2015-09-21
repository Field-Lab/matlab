% This just truns the xml-files to workable mat-files
% Chose uint8 as saving
% All 0 or 1  no  negative.. stepping away from perturbative mindset 

% This will be hacklike fo a while
% Started 2013-12-04 AKHeitman
%function prep_stimBW(string_date)
%exp_nm = string_date;

% Still assumes an alterating raster,fit raster,fit  structure 

%%% Don't have xml file writer ... but have the file names!
% requires a well maineed datarun .block.xml

%%% Backup copy of xml files exist here
%%% /netapp/snle/home/snl-e/glm/output_AH/BW_Movie/xml_interwoven
%%% on Alligator here  /Users/akheitman/xml_interwoven
function prep_stimBW(string_date)
exp_nm = string_date;

try 
    [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v19_split(exp_nm,  'BW', 'mapPRJ');
catch
    [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v19_split(exp_nm, 'BW', 'mapEI');
end

basesavedir = sprintf('%s/%s',DirPars.stimulimatxmlfiles, StimulusPars.BW.fit_rast_fullname)
% XML file names and block structures 
raster_xmlfilename = datarun_slv.block.xml{1};
SPars = StimulusPars.BW;
n_fitblock = length(SPars.NovelBlocks);
novel_xml = cell(1, n_fitblock);
for i_blk = 1:n_fitblock
    novel_xml{i_blk}.fitblock_num = i_blk;
    novel_xml{i_blk}.totalblock_num = StimulusPars.BW.NovelBlocks(i_blk);
    novel_xml{i_blk}.xmlfilename = datarun_slv.block.xml{ StimulusPars.BW.NovelBlocks(i_blk) };
end
clear i_blk

% Want the Params here to be read from the xml-files   then verify withSPARs


%% Write the Raster Movie %
clear BWmovie

mvi = load_movieAH(raster_xmlfilename, datarun_slv.block.trig{ SPars.StaticBlocks(1) });

BWmovie.moviename           = SPars.fit_rast_fullname;
BWmovie.params.height = mvi.getHeight; 
BWmovie.params.width  = mvi.getWidth;
BWmovie.params.raster_frames = SPars.raster_frames;
BWmovie.params.fit_frames    = SPars.fit_frames;
BWmovie.params.rasterblocks        = SPars.StaticBlocks;
BWmovie.params.fitblocks           = SPars.NovelBlocks;
BWmovie.params.refresh_time = mvi.getRefreshTime;
BWmovie.params.refresh_time_note = 'refreshtime is fram rate in msecs';
BWmovie.fitmovie          =[];


mov = zeros(BWmovie.params.width, BWmovie.params.height, BWmovie.params.raster_frames);
for i_frame= 1:BWmovie.params.raster_frames
    clear F
	F = mvi.getFrame(i_frame).getBuffer;
	F = reshape(F,3,BWmovie.params.width,BWmovie.params.height);   
	mov(:,:,i_frame) = F(1,:,:); 
end
% Convert to a uint binary
mov_bindouble    = round(mov);
mov_binint8      = uint8(mov_bindouble);
raster_movie = mov;

eval(sprintf('save %s/rastermovie_only.mat raster_movie', basesavedir))
% Fill in the raster movie portion
BWmovie.rastermovie.note1        = 'matrix is binary 0 1 form of the movie, saved at int8';
BWmovie.rastermovie.note2        = 'move away from perturbative framework.. no light is 0, not a negative value';
BWmovie.rastermovie.matrix       = mov_binint8;
BWmovie.rastermovie.xmlfile = raster_xmlfilename;
BWmovie.fullStimPars        = SPars;                    
clear F mvi i_frame mov mov_bindouble mov_binint8 ans

%%  Writing the fit blocks 


BWmovie.fitmovie.note1          = 'matrix is binary 0 1 form of the movie, saved at int8';
BWmovie.fitmovie.note2          = 'move away from perturbative framework.. no light is 0, not a negative value';
BWmovie.fitmovie.movie_byblock  = cell(1,n_fitblock);
BWmovie.fitmovie.ind_to_block   = SPars.NovelBlocks

for i_blk = 1:n_fitblock     
    % NEED TRIGGERS TO LOAD MOVIE  ... THERE ROLE IS NOT ENTIRELY CLEAR
    mvi = load_movieAH(novel_xml{i_blk}.xmlfilename, datarun_slv.block.trig{ novel_xml{i_blk}.totalblock_num });
	mov = zeros(BWmovie.params.width, BWmovie.params.height, BWmovie.params.fit_frames);
	for i_frame = 1:BWmovie.params.fit_frames
        F = mvi.getFrame(i_frame).getBuffer;
        F = reshape(F,3,BWmovie.params.width,BWmovie.params.height);   
        mov(:,:,i_frame) = F(1,:,:); 
    end    
    mov_bindouble    = round(mov);
    mov_binint8      = uint8(mov_bindouble);
    
    
    BWmovie.fitmovie.movie_byblock{i_blk}.xml    = novel_xml{i_blk}.xmlfilename;
    BWmovie.fitmovie.movie_byblock{i_blk}.matrix = mov_binint8;
    display(sprintf('Finished Block %d out of %d', i_blk, n_fitblock));
end

eval(sprintf('save %s/BWmovie.mat BWmovie', basesavedir))
%%
             


                 
%{                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INTERESTING NOTES ON MVI LOADING%%%
%%%    TRIGERRS ARE STRANGELY USED   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENT TRIGGER VALUES DON'T SEEM TO MATTER?
a2 =datarun_slv.block.trig{2};
a4 =datarun_slv.block.trig{4};
mdffilename = novel_xml{2}.xmlfilename;
mvi2 = load_movieAH( mdffilename , a2);
mvi4 = load_movieAH( mdffilename , a4);

mov2 = zeros(80, 40, 3600);
mov4 = zeros(80, 40, 3600);
for i = 1:3600
	F2 = mvi2.getFrame(i).getBuffer;
	F2 = reshape(F2,3,80,40);   
	mov2(:,:,i) = F2(1,:,:);
                            
    F4 = mvi4.getFrame(i).getBuffer;
	F4 = reshape(F4,3,80,40);   
	mov4(:,:,i) = F4(1,:,:); 
end                         
mdiff = mov2 - mov4;
maxdiff = max(abs(mdiff(:)))  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %%%%%%%%%%%                       
        case 'NSEM'  %% %  
            %%%%%%%%  STILL NEED TO WORK THIS OUT %%%%%%%%%
            ts_vs.width         = NSEM_Width;
            ts_vs.height        = NSEM_Height;
            ts_vs.refresh_time  = NSEM_Refresh; 
            d_mov = sprintf('%s/NSEM',output_dir);
            NovelBlock = 0;
            for k = [1,NSEM_Blocks] %%% LOAD ALL MOVIES FROM STIM_NSEM FOLDER
                if k == 1
                     load(sprintf('%s/eyemov%d',d_mov,k),'Xr'); % Xr: 3600x80x40, [0,255], double.
                     Xr = reshape(Xr,fr_sec * nsec_o ,NSEM_Width*NSEM_Height); % now 3600x3200
                     ts_vs.mov{k} = reshape(Xr',NSEM_Width,NSEM_Height,fr_sec * nsec_o );  
                     testmovie    = ts_vs.mov{k};
                end

                %%%%%%%%%% THIS SECTION IS A BIT HARD CODED DUE TO HOW THE
                %%%%%%%%%% NSEM MOVIES ARE LOADED AND SAVED INTO MATLAB!!!!
                if ~isempty ( NSEM_Blocks == k )  %%% right now using nec_e is 2x as long as nsec_o
                     display(sprintf('LOADING NSEM FIT MOVIE %d out of %d',ceil(k/2),ceil(n_blk/2)) );   
                     NovelBlock = NovelBlock + 1;
                     Xr0 = zeros(fr_sec*nsec_e,NSEM_Width*NSEM_Height);
                     load(sprintf('%s/eyemov%d',d_mov,NovelBlock),'Xr'); % Xr: 3600x80x40, [0,255], double.
                     Xr0(1:fr_sec * nsec_o,   :) = reshape(Xr,fr_sec*nsec_o,NSEM_Width*NSEM_Height); % now 3600x3200
                     NovelBlock = NovelBlock + 1;
                     load(sprintf('%s/eyemov%d',d_mov,NovelBlock),'Xr'); % Xr: 3600x80x40, [0,255], double.
                     Xr0(fr_sec * nsec_o + 1:fr_sec*nsec_e,:) = reshape(Xr,fr_sec*nsec_o,NSEM_Width*NSEM_Height); % now 2*3600 x 3200
                     ts_vs.mov{k} = reshape(Xr0',NSEM_Width,NSEM_Height,fr_sec*nsec_e); 
                end
            end

            NSEM_dir = sprintf('%s/NSEM',output_dir);
            cd(NSEM_dir);
            save('testmovie.mat','testmovie'); 
            save('tsvs.mat','ts_vs');
            cd(d_save);
    end
            %}
end