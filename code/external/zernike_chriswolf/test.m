% **********************************************************************
% Christian Wolf, http://liris.chrs.fr/christian.wolf
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with
% a note saying that the original programs are available from my web page.
%
% The programs and documents are distributed without any warranty, express 
% or implied.  As the programs were written for research purposes only, they 
% have not been tested to the degree that would be advisable in any 
% important application.  
% All use of these programs is entirely at the user's own risk.
%
% **********************************************************************

% ---- Test whether there is any difference in reconstructing with 
% ---- a full basis or a non redundant basis
% ---- There shouldn't be one, of course.

global DEB0 DEB1;
DEB0=[];
DEB1=[];

imsize=20;

load ../src-stshapes/runningman.mat
m=manrunning(11:11+imsize-1,11:11+imsize-1);
zs1=zernike_bf(imsize,4,1);
zs0=zernike_bf(imsize,4,0);
v0=zernike_mom(m,zs0);
v1=zernike_mom(m,zs1);
rec0=zernike_rec(v0,imsize,zs0);
rec1=zernike_rec(v1,imsize,zs1);
sortrows(DEB0,[ 1 2])
sortrows(DEB1,[ 1 2])
fprintf ('Difference: %f\n',sum(sum(rec0-rec1)));