function plotBundleVoltage(path, patternNos, varargin)
% PLOTBUNDLEVOLTAGE() plots the estimated voltage of an axon bundle at
% varying stimulus amplitudes
%  inputs:   path - path to the stimulus data separated by patterns
%            patternNos - patterns to plot
%            
%   optional:  display - true or false, do or don't display the array as an
%               approximate location of the bundle is found for
%               each movie
%              singleAxonVoltage - a value to be used as the 'average'
%               value of a single axon's voltage. very rough estimate
%  outputs:  
%
% usage: plotBundleVoltage('/Volumes/Analysis/2015-04-09-2/data003/', 67, 'display', false, 'singleAxonVoltage', 16.43543)

%set defaults for optionals
display = false;
singleAxonVoltage = 16.05148356;

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'display'
            display = varargin{j*2};
        case 'singleaxonvoltage'
            singleAxonVoltage = varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

bundleMeans = getBundleVoltagesAStar(path, patternNos, display);%getBundleVoltages(patternNos, display);

figure; xlabel('Stimulation amplitude (uA)'); ylabel('Average bundle voltage (mV)'); set(gcf, 'InvertHardCopy', 'off');
for patternIndex = 1:size(patternNos, 2)
    hold on; scatter(abs(bundleMeans(:, 2, patternIndex)), abs(bundleMeans(:, 1, patternIndex)));
end

ylims = ylim;
axonupper = ylims(2)/singleAxonVoltage;

legend(num2str(patternNos'), -1);

set(gca,'Box','off');   %# Turn off the box surrounding the whole axes
axesPosition = get(gca,'Position');          %# Get the current axes position
hNewAxes = axes('Position',axesPosition,'Color','none','YLim',[0 axonupper],'YAxisLocation','right','XTick',[],'Box','off');                %#   ... and no surrounding box
ylabel(hNewAxes,'Approx. # of axons');
end