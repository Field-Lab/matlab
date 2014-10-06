function javaswappath(adds, rms)
% JAVASWAPPATH  Adds some Java paths and removes others
% usage: javaswappath(adds, rms)
%
% phli 2010-05
%

javarmpath_quiet(rms);
javaaddpath(adds);