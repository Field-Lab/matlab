function [error, bad] = checkError(imm)
%CHECKERROR Plot class constructor.
%   checkError(imm) displays how many electrodes have been analyzed
%   successfully. 
%
%   tamachado@salk.edu 1/23/08

% Get list of known bad electrodes
bad = getBadElectrodes(imm.data);


% Number of electrodes analyzed
a = length(find(imm.error >= imm.SUCCESS_CLUSTER));

% Number of electrodes analzed + points filled in
b = length(find(imm.error == imm.SUCCESS_FIT));

% Number of bad electrodes
c = length(find(bad == true));

% Print information to screen
if imm.params.verbose == true
    disp('Analyzed indicates clustering has finished')
    disp('Completed indicates out of sample fitting has finished')
    disp('')
    disp(sprintf('Electrodes  analyzed: %d', a))
    disp(sprintf('Electrodes completed: %d', b))
    disp(sprintf('\nBad electrodes: %d', c))
end

% Return error array to calling function
error = imm.error;
