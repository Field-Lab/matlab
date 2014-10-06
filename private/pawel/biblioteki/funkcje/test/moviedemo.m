% VIBES  Movie of the vibrating L-shaped membrane.
%
% This is the same demo that comes with MATLAB expo, demonstrating MATLAB's
% own movie capability. Here it is also done as MPEG and gif movies. The 
% resulting movies require much less storage due to the more efficient formats.
% MPEG movies play much slower than MATLAB's own movies, GIF movies may
% require more storage in some cases.
%
% See also simplemovie

%  C. B. Moler, 6-30-91, 8-27-92.
%  Expo demo adapted by Ned Gulley, 6-21-93
%  Copyright (c) 1984-94 by The MathWorks, Inc.
%
%       Adapted to demonstrate movie package 
%       by ALW 11-27-95, 9/24/96, 2/20/02

% Demo initialization
if ~exist('MovieGUIFlag'), figNumber=0; end;

hlpStr1= ... 
        ['                                                  '
   '                                                  '
   ' Moviedemo: A demonstration of MPG and GIF        ' 
   '            movies, written from Matlab.          '
   '                                                  ' 
   '                                                  ' 
   ' This demo demonstrates the effective storage     '
   ' of the MPEG compression standard and the GIF     '
   ' movie format. Using the standard vibes example   '
   ' from the Matlab demo, we create the standard     '
   ' file (M.mat), an mpeg movie (vibes.mpg) and      '
   ' a gif movie (vibes.gif). The .mat file is        '
   ' approximately 10 times as big as the .mpg file,  '
   ' the factor is about 20 for the .gif file.        '
   ' The resolution of all movies (in pixels) is      '
   ' (almost) the same (mpeg and gif need multiples   '
   ' of 16). Avantages are selectable framesize       '
   ' and the fact that mpeg and gif movies can be     ' 
   ' played independently of Matlab.                  '
   '                                                  '
   ' press any key to continue ...                    ' 
   '                                                  ' 
   '                                                  ' 
   '                                                  '];
hlpStr2= ...      
  ['                                                  '
   ' Here some info on vibes (from the original       '
   ' Matlab demo):                                    '
   '                                                  '
         ' This movie animates a solution of the wave       '  
         ' equation for the vibrations of an L-shaped       '  
         ' membrane.  The solution is expressed as a        '  
         ' linear combination, with time-dependent          '  
         ' coefficients, of two-dimensional spatial         '  
         ' eigenfunctions. The eigenfunctions have been     '  
         ' pre-computed using the function MEMBRANE,        '  
         ' and saved in a file. The first of these          '  
         ' eigenfunctions, the fundamental mode, is the     '  
         ' MathWorks logo.                                  '  
         '                                                  '  
         ' The L-shaped geometry is of particular interest  '  
         ' mathematically because the stresses approach     '  
         ' infinity near the reentrant corner. Conventional '  
         ' finite difference and finite element methods     '  
         ' require considerable time and storage to         '  
         ' achieve reasonable accuracy. The approach        '  
         ' used here employs Bessel functions with          '  
         ' fractional order to match the corner singularity.'
         '                                                  '];
disp(hlpStr1); pause;
disp(hlpStr2); 

% This demonstration solves the wave equation for the vibrations
% of an L-shaped membrane.  The solution is expressed as a linear
% combination, with time-dependent coefficients, of two-dimensional
% spatial eigenfunctions.  The eigenfunctions have been pre-computed,
% using the function MEMBRANE, and saved in a file.  The first of
% these eigenfunctions, the fundamental mode, is the MathWorks logo.
%
% The L-shaped geometry is of particular interest mathematically because
% the stresses approach infinity near the reentrant corner.  Conventional
% finite difference and finite element methods require considerable
% time and storage to achieve reasonable accuracy.  The approach used
% here employs Bessel functions with fractional order to match the
% corner singularity.
%
% As the solution is computed and displayed at each of 12 time steps,
% snapshots of the graphics window are saved in a large matrix.
% The total storage required is proportional to the area of the graphics
% window.  About 2 megabytes is required for the default window.

% Now make the movie
% Load eigenfunction data.
load vibesdat
n = max(size(L1));
nh = fix(n/2);
x = (-nh:nh)/nh;

% Get coefficients from eigenfunctions.
clear c
for k = 1:12
     eval(['c(k) = L' num2str(k) '(24,13)/3;'])
end
 
% Set graphics parameters.
handle = figure(1);
clf
axis([-1 1 -1 1 -1 1]);
caxis(26.9*[-1.5 1]);
colormap(hot);
axis off
hold on

% Generate the movie.
delt    = 0.1;
nframes = 12;
M       = moviein(nframes,handle);
for k = 1:nframes,
    % Coefficients
    t = k*delt;
    s = c.*sin(sqrt(lambda)*t);

    % Amplitude
    L = s(1)*L1 + s(2)*L2 + s(3)*L3 + s(4)*L4 + s(5)*L5 + s(6)*L6 + ...
         s(7)*L7 + s(8)*L8 + s(9)*L9 + s(10)*L10 + s(11)*L11 + s(12)*L12;

    % Velocity
    s = s .* lambda;
    V = s(1)*L1 + s(2)*L2 + s(3)*L3 + s(4)*L4 + s(5)*L5 + s(6)*L6 + ...
      s(7)*L7 + s(8)*L8 + s(9)*L9 + s(10)*L10 + s(11)*L11 + s(12)*L12;

    % Surface plot of height, colored by velocity.
    % Cut out the reentrant corner
    V(1:nh,1:nh) = NaN*ones(nh,nh);
    cla
    surf(x,x,L,V);
    axis off
    
    % ALW 11-27-95, 9/24/96
    % Save a frame
    drawnow
    M(:,k) = getframe(handle);
    disp(['Saving frame ',num2str(k)])
    makeframe(k,'vibes',[560,400]);     
    % Size of Standard Matlab window is(was) "golden ratio" 560,396 
    % but we need multiples of 16 for the mpeg encoder.
    % END ALW
end;
hold off
%====================================
% Store the movie
mvstore(figNumber,M);

% ALW 11-27-95, 9/24/96, 05-15-98, 02-20-02
makemovie([1:nframes],'vibes')
save M M

disp('Now playing 10 repetitions of matlab movie (under Matlab)')
clf
axis off
movie(M,10)

pdir = makemovie;
disp(' ')
disp('Now playing endless loop of mpeg movie.')
disp('Kill window (Click on Exit if applicable) to continue with moviedemo.')
unix(['xterm -e ',pdir,computer,'/mpeg_play -loop ',pwd,'/vibes.mpg']);

disp(' ')
disp('Now playing gif movie: requires netscape on your system.')
disp('Demo may be slow due to image size, see simplemovie for a faster demo.')
disp(' ')
disp('Exit netscape to continue with moviedemo.')
unix(['netscape file://',pwd,'/vibes.gif']);


disp(' ')
disp('===================')
disp('done with moviedemo')
disp('===================')
% END ALW
