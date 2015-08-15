function [pats_yes, pats_no] = checkForElecResp(fpath, neuronId, pattern)
% Checks to see if a given elecResp file exists. 
% input: fpath 
dirInfo = dir(fpath);
for nn = 1:length(neuronId);
    neuron = neuronId(nn);
    for  n = 3:size(dirInfo,1) %%106 %
        if ~dirInfo(n).isdir && ~strcmp(dirInfo(n).name,'activationResults.mat')
            fname = dirInfo(n).name;
            i = find(fname=='_',2,'first');
            if strcmp(['n' num2str(neuron)],fname(i(1)+1:i(2)-1))
                temp = load([fPath fname]);
                elecResp = temp.elecResp; clear temp;
            end
        end
    end
end
end %function end