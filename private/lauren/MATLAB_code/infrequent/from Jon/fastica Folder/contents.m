% FastICA for Matlab 5.x
% Version 1.01  April 23 1998
% Copyright (c) Jarmo Hurri, Hugo Gävert, Jaakko Särelä, and Aapo Hyvärinen.
%
% Type fasticag to launch the graphical user interface
%
%
% FastICA programs:
%   fasticag  - Graphical user interface for FastICA
%   fastica   - command line version of FastICA
%
% Separate functions used by FastICA programs.
%   fpica     - main algorithm for calculating ICA
%   whitenv   - function for whitening data
%   pcamat    - calculates the PCA for data
%   dispsig   - plots the data vectors
%   remmean   - function for removing mean
%   selcol    - function used by pcamat
%
%   gui_cb    - needed by fasticag
%   gui_adv   - needed by fasticag
%   gui_advc  - needed by fasticag
%   gui_l     - needed by fasticag
%   gui_lc    - needed by fasticag
%   gui_s     - needed by fasticag
%   gui_sc    - needed by fasticag
%   gui_cg    - needed by fasticag
%   gui_help  - needed by fasticag
%
% Misc.
%   demosig   - generates some test signals
