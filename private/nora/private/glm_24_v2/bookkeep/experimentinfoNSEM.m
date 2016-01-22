function [exp_info] = experimentinfoNSEM(string_piece) 
% MY_FUNCTION   Holds onto minimal information of datarun and stimulus type
%               for the natural scene and decomposition pieces
%
% arguements:   string_piece - string of experiment piece with NSEM runs
%                                ie. '2012-08-09-3'    
%            
% output:       exp_info - full sructure with minimal information
%                   .dr  - name of each datarun  ie data000
%
%
% AKHeitman 2014-04-02


%%% IMPORTANT %%%
% Needs to be updated for each new NSEM piece added

%%% NOTES %%%
% NSEM  : Natural Scenes with Eye Movements  fit and raster
% Decomp: Decomposed Eye Movements
% WN    : BW White noise intrleaved fit and raster
% mas   : master, usually an cell-typing RGb run
% slv   : slave of the master

%%% UPDATES %%%


switch string_piece
    case '2012-08-09-3'
        exp_info.dr.mas          = 'data002';
        exp_info.dr.slvNSEM      = 'data005';
        exp_info.dr.slvWN        = 'data006';
        exp_info.dr.slvDecomp    = 'data004';
          
        exp_info.masfile         = 'RGB-8-1-0.48-11111-80x40';
        exp_info.WNfile          = 'BW-8-1-0.48-11111_RNG_16807';
        exp_info.NSEMfile        = 'NSEM_eye-120-3_0-3600'; 
        exp_info.NSEMscheme     = 'schemeA';
        exp_info.Decompfile      = '';
        exp_info.DecompScheme    = '';
        
    case '2012-09-27-3'
        exp_info.dr.mas      = 'data003';
        exp_info.dr.slvNSEM  = 'data002';
        exp_info.dr.slvWN    = 'data005';
        exp_info.dr.slvDecomp= 'data004';
      
        exp_info.masfile         = 'RGB-10-2-0.48-11111-64x32';
        exp_info.WNfile          = 'BW-8-1-0.48-11111_RNG_16807';
        exp_info.NSEMfile        = 'NSEM_eye-120-3_0-3600'; 
        exp_info.NSEMscheme     = 'schemeA';
        exp_info.Decompfile      = '';
        exp_info.DecompScheme    = '';
        
    case '2013-08-19-6'
        exp_info.dr.mas      = 'data000';
        exp_info.dr.slvNSEM  = 'data001';
        exp_info.dr.slvWN    = 'data003';
        exp_info.dr.slvDecomp= 'data002';
        
        exp_info.masfile         = 'RGB-10-2-0.48-11111-64x32';
        exp_info.WNfile          = 'BW-8-1-0.48-11111_RNG_16807';
        exp_info.NSEMfile        = 'NSEM_eye-long-v2'; 
        exp_info.NSEMscheme     = 'schemeA';
        exp_info.Decompfile      = 'NSEM_eye-trace-2';
        exp_info.DecompScheme    = '';
        
    case '2013-10-10-0'
        exp_info.dr.mas      = 'data000';
        exp_info.dr.slvWN    = 'data005';
        exp_info.dr.slvNSEM  = 'data001';
        exp_info.dr.slvDecomp= 'data002';
        
        exp_info.masfile         = 'RGB-10-2-0.48-11111-64x32';
        exp_info.WNfile          = 'BW-8-1-0.48-11111_RNG_16807';
        exp_info.NSEMfile        = 'NSEM_FEM900FF_longrast'; 
        exp_info.NSEMscheme     = 'schemeB';
        exp_info.Decompfile      = 'NSEM_decomp_5TR60IM';
        exp_info.DecompScheme    = '';
        
     case '2014-06-04-0'
        exp_info.dr.mas      = 'data004';
        exp_info.dr.slvWN    = 'data003';
        exp_info.dr.slvNSEM  = 'data002';
        exp_info.dr.slvDecomp= '';
        
        exp_info.masfile         = 'RGB-10-2-0.48-11111-64x32';
        exp_info.WNfile          = 'BW-8-1-0.48-11111_RNG_16807';
        exp_info.NSEMfile        = 'NSEM_FEM900FF_longrast'; 
        exp_info.NSEMscheme     = 'schemeB';
        exp_info.Decompfile      = 'NSEM_decomp_5TR60IM';
        exp_info.DecompScheme    = '';
end


%{
stim.master.pixelsize = 10;
stim.master.height = 32; stim.master.width  = 64;
stim.master.refreshrate = 2;
stim.master.refreshnote = 'rfreshrate is units of (1/120) seconds between stim refresh';
stim.master.type = 'RGB';
if strcmp(exp_nm , '2012-08-21-1') || strcmp(exp_nm , '2012-08-09-3') || strcmp(exp_nm, '2012-04-13-4')
    stim.master.height = 40; stim.master.width  = 80; stim.master.pixelsize = 8;
end

if strcmp(exp_nm, '2012-08-09-3')
   stim.master.refreshrate = 1;
end
%}

%{
if strcmp(string_date, '2013-10-10-0'); dn_mas = 'data000';             dn_BW = 'data005';   dn_NSEM = 'data001'; end
if strcmp(string_date, '2013-08-19-6'), dn_mas = 'data000');    dn_BW = sprintf('data003'); dn_NSEM = sprintf('data001'); end
if strcmp(string_date,'2012-09-27-3'),  dn_mas = 'data003');    dn_BW = sprintf('data005'); dn_NSEM = sprintf('data002');  end
if strcmp(string_date,'2012-08-09-3'),  dn_mas = 'data002');    dn_BW = sprintf('data006'); dn_NSEM = sprintf('data005');  end
if strcmp(string_date,'2012-09-21-1'),  dn_mas = 'data004');    dn_BW = sprintf('data006'); dn_NSEM = sprintf('data005');  end
if strcmp(string_date,'2012-09-24-2'),  dn_mas = 'data002');    dn_BW = sprintf('data005'); dn_NSEM = sprintf('data001');  end
if strcmp(string_date,'2012-09-13-1'),  dn_mas = sprintf('data004');    dn_BW = sprintf('data007'); dn_NSEM = sprintf('data003');  end
if strcmp(string_date,'2012-04-13-4'),  dn_mas = sprintf('data000');    dn_BW = sprintf('data003'); dn_NSEM = sprintf('data002');  end
if strcmp(string_date,'2012-01-27-4'),  dn_mas = sprintf('data007');    dn_BW = sprintf('data003'); dn_NSEM = sprintf('data005');  end
%}

%{
schemeA.rasterblocks = 1:2:119;
schemeA.fitblocks   = 2:2:120;
schemeA.rasterblockframes = 3600;
schemeA.rasterblocksecs   = 30;
schemeA.rasterblocktriggers    = 37;
schemeA.fitblockframes = 7200;
schemeA.fitblocksecs   = 60;
schemeA.fitblocktriggers  = 73;
schemeB.rasterblocks = 1:2:59;
schemeB.fitblocks   = 2:2:60;
schemeB.rasterblockframes = 14400;
schemeB.rasterblocksecs   = 120;
schemeB.rasterblocktriggers = 145;
schemeB.fitblockframes = 14400;
schemeB.fitblocksecs   = 120;
schemeB.fitblocktriggers  = 145;
schemeB.note = 'roughly 5 second reset between blocks';
if strcmp(dn.NSEMscheme, 'schemeA')
    fitscheme = schemeA;
end
if strcmp(dn.NSEmscheme, 'schemeB')
    fitscheme = schemeB;
end
%}


%{
%%%%% BW SLV  %%%
stim.BW.commonmovie         = true;
if strcmp(exp_nm, '9999-99-99-9')
    stim.BW.commonmovie = false;
end
stim.BW.default_tstim = .00832750;
stim.BW.pixelsize = 8; %pixelsize 
stim.BW.height = 40;
stim.BW.width = 80;
stim.BW.refreshrate = 1; %refresh rate
stim.BW.framespertrigger = 100; % frame per trigger
stim.BW.seedS = 11111; % stat seed 
stim.BW.seedA = 16807; %random number generating component A
stim.BW.seedM = 2^31 -1; % random number generating component M
stim.BW.seedC = 0  ; % rando
stim.BW.fit_frames = 3600;
stim.BW.raster_frames = 1200;

if strcmp(exp_nm, '2012-01-27-4')
    stim.BW.seedA = 1664525;
    stim.BW.seedC = 1013904223;
    stim.BW.seedM = 2^32;
    stim.BW.seedS = 11111; % initial seed
    stim.BW.blk = 2*[11,16,17,18,21,22,23,24,25,26,27,30,31,34,35,36,37,38,41,43,45,46,47,48,51,53,54,55,57,58];
end
        
  %%% 2013-10-10-0  NSEM - stimulu
  ;data001
;trigger 125s duration: 7800
(let* ((display *display*)
       (frames 14400)
       (repeats 30) ;; number of trials
       (path "Copperdome:Stimuli:eyemovement:FEM900FF_longrast.rawMovie")
       (wrap (- (* 3600 120) frames)) ;; wrap around before hitting this length of movie
       (common (list :type :movie :path path :buffer-size 100 :interval 1 :back-rgb #(0.23 0.23 0.23)
                     :stixel-width 2 :stixel-height 2 :x-start 0 
                     :x-end 640 :y-start 80 :y-end 400))
       (frames-A 14400)
       (frames-B (* frames-A 1))
       (start-frame-A 0)
       (start-frame-B frames)
       (stimulus-A (apply 'make-stimulus :start-frame start-frame-A :frames frames-A common))
       (stimulus-B (apply 'make-stimulus :start-frame start-frame-B :frames frames-B common))
       (wait-trigger t)
       )
  (prepare-display display stimulus-A) ;; spatial pattern is unchanged, so prepare both for speed
  (prepare-display display stimulus-B)
  (without-interrupts
   (dotimes (i repeats)
     (when (plusp (abort-button-down)) (abort))
     (animate-display display stimulus-A :wait-trigger wait-trigger :beep nil)
     (animate-display display stimulus-B :wait-trigger wait-trigger :beep nil)
     (incf (start-frame stimulus-B) frames-B)
     (when (>= (start-frame stimulus-B) wrap)
       (setf (start-frame stimulus-B) frames))
     ))
  )
%}

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTEBOOKS BELOW  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%{
2011-12-04-0   33 degrees
AH Comments
-looks like potentially bad animal
Notebook
data000 RGB 10-16-0.48-11111
data001 DS sin (didn't stop recording befor seting up new stimulus - has some fake trigger at the end)
data002 DS dots  black
data003 DS dots  black and white
GLM eyemovent/white noise/jitter
data004 BW 8-1-.0.48- 30s repaets frozen and movig seed, 50 reps 
data005 eye movments repaets frozen and movig seed, 50 reps 
data006 BW 8-1-.0.48 jitter true 30s repaets frozen and movig seed, 50 reps 
data007 RGB 10-16-0.48-11111
data006 BW 16-1-.0.48 true 30s repaets frozen and movig seed, 30 reps
data009 eye movments shift by50px repaets frozen and movig seed, 50 reps
data010 RGB 10-16-0.48-11111

2011-12-13-1  35 degrees
AH Comments
Notebook
data000 RGB 8-2-0.48-11111
data001 RGB 4-2-0.48-11111
data002 DS moving sinusoid eight direction
data003 DS moving dots 
data004 GLM BW 8
data005 GLM eye-movement
data006 GLM BW 8 jitter
data007 RGB 8-16-0.48-11111: streaming is possibly wrong
data008 saccade movie (gray->saccade->natim): stimlus frozen in the middle
data009  RGB 8-16-0.48-11111: found out that the retina should be degraded in the last couple of runs (starting data007; could be earlier)
Martin did imaging for Peter.  w/ Penut Lectin.

2012-01-27-0  (35 degrees)
AH Comments
Notebook
data003 eye repeats 5 images, 2 fliped, 4 traces, 1.5s, 50rep, 3000s  STOPPED.
note: array was moved at this point.
data004 eye repeats 5 images, 2 fliped, 4 traces, 1.5s, 50rep, 3000s  RERUN. -- around 2270 the stimulus machine frozen -- TERMINATED.
path "Chinacat:Stimuli:eyemovement:eye-traces.rawMovie")
-has NSEM fitting data as well
-not sure if this is the same eye-trace.rawMovie that we have now
%%%%
data000 RGB 8-2-0.48-11111
data001 DS sin
data002 DS dots
data003 eye repeats 5 images, 2 fliped, 4 traces, 1.5s, 50rep, 3000s  STOPPED
note: array was moved at this point.
data004 eye repeats 5 images, 2 fliped, 4 traces, 1.5s, 50rep, 3000s  RERUN. -- around 2270 the stimulus machine frozen -- TERMINATED.
recreate stinulus file last entries
data005 RGB 8-2-0.48-11111, 900s.
data006 BW repeats -- BW -8-1-.48-11111 block design, GLM fitting, 60 rep, frozen odd block 10, random seeds, 30 
data007 BW repeats -- BW -8-1-.48-11111-jitter block design, GLM fitting, 60 rep, frozen odd block 10, random seeds, 30 
data008 NSEM repeats -- aka eye movement data, block design, GLM fitting, 60 rep, frozen odd block 30, random seeds, 60 
data009 NSEM repeats -- restarted the above in 5 sec -- 3000+ sec, and stop
data010 RGB 8-2-0.48-11111, 900s. 

2012-01-27-4 (35 degrees)
AH Comments
Notebook
data000 RGB 8-2-0.48-11111
data001 DS sin
data002 eye repeats 5 images, 2 fliped, 4 traces, 1.5s, 50rep, 3000s, stopped 3100s.
data003 BW repeats -- BW -8-1-.48-11111 block design, GLM fitting, 60 rep, frozen odd block 10, random seeds, 30 
data004 BW repeats -- BW -8-1-.48-11111-jitter block design, GLM fitting, 60 rep, frozen odd block 10, random seeds, 30 
data005 NSEM repeats -- aka eye movement data, block design, GLM fitting, 120 rep, frozen odd block 30, random seeds, 60, stopped around 8500 sec (stimulus machine frozen)
data006 RGB 8-2-0.48-11111, 1800s. -- file not readable
data007 RGB 8-2-0.48-11111, 1800s. -- rerun
 
2012-04-13-0  (37 degrees) 
AH Comments
Notebook
dat000: NSEM, 10 images, 100 repeats (~3pm)
dat001: Classification (~3:20pm -- 3:54pm)      RGB 8-2-0.48-11111  (80 by 40)  
dat002: NSEM block-design for GLM (a total of 90min; 3:57pm --5:32pm): in the lisp code, the number of repeats was 120 instead of 60: we stopped 5673 elapsed time.
dat003: BW block-design 8x8 (a toal of 40min; 5:35pm -- 6:18pm)
dat004: BW jitter block-design 8x8 (a total of 40min; 6:19pm -- 6:59pm)
dat005: Classification (30min; 7:00pm -- 7:30pm)
dat006: NSEM de-constructed stimulus (7:35pm -- ~9:10pm; elapsed time: 5712 sec)

2012-04-13-4 (37 degrees)
AH Comments
nice On and Off Parasol mosaics  even SBCs
60 micron board ..  this looks promising
will proceed to analyze the repeats

Need to start moving over the NSEM runs !! ! ! 
Notebook
data000: RGB 8-2-.48-11111 Classification (30min; 11:15PM -- )
data001: NSEM, 10 images, 100 repeats (17.5min; 11:47PM -- )
data002: NSEM block-design for GLM (a total of 90min = 5490 sec; ~12:10am -- 1:46am) elapsed time: 5562 sec
data003: BW block-design 8x8 (a toal of 2424 sec; 1:49am -- ~2:31am)
data004: BW jitter block-design 8x8 (a total of 2424 sec; 2:32am -- 3:13am)
data005: Classification (30min = 1800sec; 3:14am -- 3:44am)
data006: NSEM de-constructed stimulus (max 2 hrs: 3:47am -- 
data007: NSEM, 10 images, 100 repeats (17.5min = 1050sec; 6:05AM -- 6:22AM)
data008: Classification (30min; 6:23AM -- stopped at1433 sec elapsed time, due to the lack of disk space)

2012-08-09-3 (38 degrees starting at data002)
AH Comments
30 micron piece

Notebook
First attempt at suction without plug.  Did not work; filter paper for drying ended up sticking to piece.  Plug put back on.
Squashed very slightly on dissection stage, then moved to rig and squashed more with cracker.  Good coverage but very small spikes.  Parasol mosaic in online looks pretty central.
data000 BW-20-5-0.48-11111-30x30-60.35 NDF 0.6
data001 BW-3-5-0.48-11111-200x200-60.35 NDF 0.6
data002   RGB-8-8-0.48-11111
MG took over from EJ and PL.  Switched from OLED to CRT.  After focusing, temperature slowly raised from 35 to 38.  Once at 38, immediately started the RGB run.  Spikes were weak (presumably due to centrality of the piece).  However, the coverage was good and the recording was extremely stable.  Good mosaics.
data003 NSEM repeated runs of movie 
10 second movie repeats of Natural Scenes with microsaccades,  10 images, one second each.   11 second trigger.  100 repeats.
data004 NSEM eye traces
-1 hour of move   I can't remember right now.
data005 NSEM 
Hour and a half of interwoven NSEM fitting, and NSEM testing.  
data006    8 by 8
Coarse BW noise run.  40 minutes.   We thought it would last longer so we recorded for close to 50 minutes... last 10 minutes was just a blank gray screen.
data007  2 by 2
Fine BW noise run.  Should go for 66 minutes


2012-08-21-1 (38 degrees)
AH Comments
(the piece me and DA co-ran)
- not super steable.. ok mosaics.. def. a hole
- but all 4 populations
- 30 micron board
-already did map analysis
- DONE MECHANICS 

Notebook
switch to CRT
turret in with half-silvered mirror (normal for CRT)
raise temp to 38
no movement of stage
data003  Bookend RGB 8-1-0.48-1111-80x40   1800 s
-very stable
data004  NSEM 60 reps of 30 sec frozen, 60 sec novel 
-upper gravity drop slowly sank from 460 to 440
data005 BW 8-1 noise, 60 reps of 10 secs frozen, 30 secs novel
data 006   BW jitter 8-1, 100 reps of 10 secs frozen, 30 secs novel 
data007 Eye Movement-Tracers  60 blocks of 1.5 sec
data008 Bookend   RGB 8-1-0.48-1111-80x40   1800 s



2012-09-13-1 (37.5 degrees .. killed the ONs)
AH Comments
Notebook
data000 RGB 10-2-.48-11111
data001 map from 000 lower corner, 3600s BW,intervall1 
Tried to raise temperature to 37.5   But the coverage and spike rate declined significantly in the process.  Brought temperature back down to 36.   Retina stabilized but with notably worse coverage and spike rate as compared to data000 and data 001.
data002 RGB 10-2-.48-11111
Streamed data, found that we have lost ALL On Parasols and and ON midgets.  Still have strong Off-parasol mosaic and will proceed with Natural scenes as well as tracers.
data003  NSEM   60 repeats of 30s frozen, 30 sec novel
-was stable thourghout the nearly 2 hour stimulus 
data004 RGB 10-2-.48-11111
data005  Eye traces
data006  RGB 10-2-.48-11111   900 sec
data007  BW 8- 1 .48 - 11111  60 repeats but finished 5 minutes early so maybe only 50-55 repeats
-10 sec frozen , then 30 sec novel

2012-09-21-1 (35 degrees)
AH Comments
-moving over NSEM runs and doing Map Analysis
- basically will ignore this
-checking white noise runs
-fairly stable!  not flat.. but ok.. 
-good on mosaics with slight hole in center  (On Midget and On Parasol!)
- offs ren't that good
Notebook
data000 Binary RGB 10-16- 0.48 - 11111
data001 NSEM Repeats (60 reps of 30 sec fixed   60 sec fi)
stoped recording error - continiue stimulus in new run
data002 NSEM Repeats (60 reps of 30 sec fixed   60 sec fi)
data003 NSEM traces
-piece was fairly stabe
data004Binary RGB 10-16- 0.48 - 11111
-still had nice good full on mosaics
data005NSEM Repeats
-full run ,  no hitch,  very slow moderate decline in firing rate
data006BW repeats
-looked good
data007BW jitter repeats
-before this was started the tissue was not reponsive to 2 by 2, weak response to 8 by 8
-still robust light response
data008 RGB bookened
-already starting to oscillate slightly

2012-09-21-3 (37.5 degrees)
AH Comments
-moving over NSEM runs and doing Map Analysis
-checking white noise runs
-rather unstable during white noise runs?. but 4 solid mosaics.. not perfect but pretty good  (for data004)
-STA circles the outlines poorly.. actually better than it looks
-60 micron board
-checking white noise runs

-lets move over the NSEM traces and repeats
Notebook
data000 Binary RGB 10-16- 0.48 - 11111
37.5C
data001 Binary RGB 10-16- 0.48 - 11111
data002 NSEM Repeats (60 reps of 30 sec fixed   60 sec fi)
data003 NSEM traces
data004Binary RGB 10-16- 0.48 - 11111
-first 600 was unstable ... then it really flattened out nicely
data005  BW 8 by 8 run
-first 500 oscilate then stable from 500 till 1700.  Then some instabilit.  Presumably due to low solution levels in the top chamber.  We are not sure how it dropped so much.  Luckly it doesn't seem to have gone dry.  Just down to roughly to 200 ml.
data006 NSEM
We seemed to have more midgets so we ran the NSEM again.  This time only 60 repeats rather than near 100.
data007 RGB 8-16 ...  30 min
data008 Lauren's stim200 


2012-09-24-2 (35 degrees and 37.5)
AH Comments  lost notebook reconstructed
- Don't worry about this yet
- 60 micron board
-OK mosaics..have almost 4 populations of ~70% mosaics.
-should be able 
-not particularly stable either
-not figure worthy.. but should check mosaics
Notebook
35C 
4.5x objective
CRT
lost notebook reconstructed
data000  RGB-10-16-0.48-11111
data001 NSEM Repeats (60 reps of 30 sec fixed   60 sec fi)
37.5C 
data002  RGB-10-16-0.48-11111
data003 NSEM traces
data004  RGB-10-16-0.48-11111
data005 BW repeats
data006 BW jitter repeats
data007  RGB-10-16-0.48-11111



2012-09-27-3  (NSEM at 37.5 ) 
AH Comments
Notebook
VERY GOOD NSEM PIECE WITH GOOD MOSAICS
data000 Binary RGB 10-2-0.48-1111-64x32
-really wierd electrical noise and some general instability
-not sure that we have hammered it out.. but the streamed data shows nice mosaics
piece suffered from poor plug lowering and prolly not squashed enough.. there is a slight gradient.. despite that the datarun looks good in visoin and we will heat up
-Gradient
Heating up tp 37.5
data001 NSEM
-got tripped up with some lab view error...  steady increase in spike rate
data002  NSEM 60 repeats
-the coverage got steadily better 
-gradient decined significantly.. piece looks much better in terms of even coverage
data003  Binary RGB 10-2-0.48-1111-64x32
-still retain the on and off parasol mosaic
data004 eye traces
-very stable
data005   BW 8-1
-at 900 secs started to tail off in spiking for no apparent reason, but seems to have stabilized by 1100
-piece slowly began to oscillate
data006  BW 32-1
-fairly stable
-same as BW 8-1  but just four times as large on each dimension      test correlated input
data007    RGB 10-2-0.48-1111-64x32   30min
data008 RGB 10-2-0.48-1111-64x32
start wash in 1030
then carbanoxenol 25uM
ran out of gas , changed the bottle
1030-2150
for short duration many noise events.. perhaps too many bubbles.  Pressure was tpp high, then reduced.
(2000 - 2500).  This seemed to fix the noise events. 
5100-5300 emty trap, liky that insability over last 1000s was inflow touchnig water level 
data009 RGB 10-2-0.48-1111-64x48



2013-06-18-0  (mostly done at 35 degrees)
datarun000: RGB-10-2-.48-11111-64X32
-not perfectly stable  but piece looks good
-very perihperal.. probaby explains the low spike rate.

datarun001: decomposition (AH version  4 sets of 10 seconds randomly ordered ... 60 times )
noise event did occur
otherwise is running smoothly 
top container got up to 650 for a little bit 


datarun002: tracers (MG version
datarun003  RGB-10-2-.48-11111
just for a quick diagnostic
slight hickup with lab view
maybe a gradient developing
just not that stable
datarun004 NSEM fitting for GLM !
datarun005:  BW GLM fit
- using size 10 instead of 8 since the receptive fields are so large
testing temperature
Temperature brought up to 38!  (sowly over an hour.. no recording)
datarun006:  RGB 10-2-.48-11111
-see response at 38 degrees

2013-06-18-7 (38 degreees)
Already 38
datarun000  RGB 10-2-0.48-11111
fairly stable.. not perfect
datarun001  Decomp  10secs at a time AH version
datarun002  ECpom  Martin's orginal 1 sec at a time verions
datarun003 NSEM GLM fit + raster
datarun004 BW GLM fit and raster
~maybe an hour gap here with AH sleeping 
datarun005 RGB 10-2-.48-11111  just for 900 secs to get diagnostics
-spike rate considerably lower
-something wierd about the ames on top.. will clean out.

2013-08-19-0 (NSEM at 38 ) 
AH Comments
Notebook
datarun000: RGB-10-2-.48-11111-64X32 35.7 degrees
-not perfectly stable  but still running
datarun001: RGB-10-2-.48-11111-64X32 38 degrees
-do map analysis from here? th real run! 
-fairly stable  .. not great
datarun002: tracers .. new version 2!!  eye-traces-2.rawMovie  38 degrees

;;; data002 ;;;;;
;; trigger 105 seconds 
(let* ((display *display*)
       (frames 12000) ;; check against movie frames
       (repeats 50)
       (stimulus (make-stimulus :type :movie
                                :path "Copperdome:Stimuli:eyemovement:eye-traces-2.rawMovie"
                                :back-rgb #(0.23 0.23 0.23)
                                :buffer-size 80 :interval 1
				:x-start 0 :x-end 640 :y-start 80 :y-end 400
                                :stixel-width 2 :stixel-height 2 :frames frames
                                )
                 ))
  (prepare-display display stimulus) ; do this only once to save time
  
  (without-interrupts
   (dotimes (i repeats)
     (when (plusp (abort-button-down)) (abort))
     (animate-display display stimulus :wait-trigger t :wait-key nil)))
  )

datarun003: NSEM fit data   old version  eye-120-3_0-3600.rawMovie  38 degrees
datarun004: BW GLM fit     38 degrees
datarun005 RGB-10-2-.48-11111-64X32

2013-08-19-3 (NSEM at 37.5 )    will work on if necessary     no work done yet
AH Comments
work on later.. since the stiching together will be a bitch !! !! !!
Notebook
datarun000  RGB 10-2-0.48-11111  37.4
datarun001 eye-long-v2.rawMovie
-new NSEm FIT !!!!
-looks pretty good!
datarn002.  eye-traces-2.rawMovie
-some labview errorkilled it early.. looked pretty good otherwise
datarun003   eye-traces-v2.rawMovie
-labview died  bt atleat i got about a half hour
datarun004 BW glm fit
-lab view died 20 minutes in .. wtf?
datarun005 RGB-10-2-0.48-111111
bookend
not stable now

%%% I THINK THE saved stimulus might get mixed up here !! !! !!

2013-08-19-6 (NSEM at 35.5 ) 
AH Comments
Notebook
datarun000  RGB-10-2 -0.48-11111  (35.5)
datarun001 eye-long-v2.rawMovie
excellent stability!
It seems to have a high tonic firing rate though... so not as good as originally thought perhaps
datarun002  eye-trace-2.rawMovie
datarun003 BW GLM fit
smooth and stable
datarun004 bookend rgb 10-2-0.48-11111
(will heat up after first 15 min)
ater 15 minutes.. slowly raised to 37.5
datarun005 eye-trace-2.rawMovie   37.5 degrees
-comparing 35.5  to 37.5




2013-10-10-0    look should be good piece!  35 degrees
AH Comments
-this looks like a potentially good piece. Full on -parasol and on midget mosaic .. almost full off parasol.  
-extremely stable piece 
-
Notebook
data000: RGB 10-2-0.48-11111 
data001: FEM900FF_longrast.rawMovie
data002  decomp_5TR60IM.rawMovie
data003:  quick RGB run 10-2-0.48-11111 
BWglmfit  data005
data006 cook RGB 10-2-0.48-11111

2013-10-10-03    36.5 degrees
AH Comments
-Full on parasol and on midget .. almost full Off Parasol
-Very large receptive field sizes..   probably a little too peripheral 
-slightly unstable .. but not horrible
Notebook
data000: RGB 10-2-0.48-11111 
data001: FEM900FF_longrast.rawMovie
data002  decomp_5TR60IM.rawMovie
data003:  quick RGB run 10-2-0.48-11111 
done

2013-10-10-6  35 degrees
AH Comments
-decent on and off parasol mosaics.. near full SBC mosaic 
-not super sttable.. firing rate slowly took off ..  but still a good piece 
Notebook
data000: RGB 10-2-0.48-11111 
data001: FEM900FF_longrast.rawMovie
data002  decomp_5TR60IM.rawMovie
data 003 Bookend  just 900 secs   
rgb-10-2-048-11111



2013-10-15-0    35 degrees
AH Comments
-this looks like a potentially good piece. Full on -parasol and on midget mosaic .. almost full off parasol.  
-extremely stable piece 
-
Notebook
data000: RGB 10-2-0.48-11111 
data001: FEM900FF_longrast.rawMovie
data002  decomp_5TR60IM.rawMovie
data003 RGB-10-2-0.48-11111

2013-10-15-2    35 degrees
AH Comments
-this looks like a potentially good piece. Full on -parasol and on midget mosaic .. almost full off parasol.  
-extremely stable piece 
-
Notebook
data000: RGB 10-2-0.48-11111 
data001: FEM900FF_longrast.rawMovie
data002  BW GLM fit  
data003 RGB-10-2-0.48-11111


2013-10-15-4    35 degrees
AH Comments
-this looks like a potentially good piece. Full on -parasol and on midget mosaic .. almost full off parasol.  
-extremely stable piece 
-
Notebook
data000: RGB 10-2-0.48-11111 
data001: FEM900FF_longrast.rawMovie
data002 RGB 900 seeconds
data003 decomp_5TR60IM.rawMovie
%}


end