function A=mix_matrix
% MIX_MATRIX   Grabs the mixing maxtrix from the screen for ICADEMO
%

% Grab GUI
A_11 = get(findobj('Tag','X_11'),'String');
A_12 = get(findobj('Tag','X_12'),'String');
A_21 = get(findobj('Tag','X_21'),'String');
A_22 = get(findobj('Tag','X_22'),'String');

% Construct mixing matrix
A = [str2num(A_11) str2num(A_12); str2num(A_21) str2num(A_22)];
