function allThresholds = plotAmplitudeModulation(modulationStatsFolder, outputFolderFigures, varargin)
% plotAmplitudeModulation(modulationStatsFolder, outputFolderFigures)
%
% Parameters:
%   - modulationStatsFolder
%   - outputFolderFigures

%% Read optional inputs, format the imput folders

imageFormat = 'epsc';
neuronList = [];

% Checking the optional parameters
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn', ...
        'Unexpected number of arguments');
    throw(err);
end

% Reading the optional input arguments
for kk=1:(nbin/2)
    if ~ischar(varargin{kk*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{kk*2-1})
        case 'imageformat'
            imageFormat = varargin{kk*2};
        case 'neuronlist'
            neuronList = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% If output folder doesn't exist, make them here
if ~exist(outputFolderFigures,'dir')
    mkdir(outputFolderFigures);
end


%% Finding all the neuron data

contentsModulationStatsFolder = dir(modulationStatsFolder);
neuronNames = struct('name','');
nNeurons = 0;
% 
% for kk=1:length(contentsModulationStatsFolder)
%     if strfind(contentsModulationStatsFolder(kk).name,'.mat')
%         nNeurons = nNeurons + 1;
%         neuronNames(nNeurons).name = contentsModulationStatsFolder(kk).name;
%     end
% end

for kk=1:length(contentsModulationStatsFolder)
    if strfind(contentsModulationStatsFolder(kk).name,'.mat')
        if ~isempty(neuronList)
            cNeuronID = str2double(contentsModulationStatsFolder(kk).name(7:end-8));
            if nnz(neuronList==cNeuronID)>0
                nNeurons = nNeurons + 1;
                neuronNames(nNeurons).name = contentsModulationStatsFolder(kk).name;
            end
        else
            nNeurons = nNeurons + 1;
            neuronNames(nNeurons).name = contentsModulationStatsFolder(kk).name;
        end
    end
end

%% For each of those neurons, plot the amplitude of the modulation

fh = figure();
hold on

allThresholds = zeros(nNeurons, 2);
for kk=1:nNeurons
    % Load the modulation data
    load(fullfile(modulationStatsFolder, neuronNames(kk).name));
    cNeuronID = str2double(neuronNames(kk).name(7:end-8));
    
    % Plot
    allThresholds(kk,1) = cNeuronID;
    allThresholds(kk,2) = plotSingleNeuronModulation(alternationData, fh, cNeuronID);
    
    saveas(fh,fullfile(outputFolderFigures, neuronNames(kk).name(1:end-8)),imageFormat);
end

end % plotAmplitudeModulation

function [threshold] = plotSingleNeuronModulation(alternationData, fh, neuronID)
% Plots the neural response to contrast reversal for gratings of increasing
% size. Also fits a function to the response.

x = alternationData.gratings;
if size(x,2)>size(x,1)
    x = x.';
end
y = alternationData.response;
if size(y,2)>size(y,1)
    y = y.';
end
e = alternationData.variances;
if size(e,2)>size(e,1)
    e = e.';
end

figure(fh); clf; set(fh,'color','white');
hold on
errorbar(x,y,sqrt(e),...
            'MarkerFaceColor',[0 0 0],'Marker','square',...
            'LineStyle','none','Color',[0 0 0]);
grid on
ylabel('#spikes elicited/reversal','fontsize',12)
xlabel('grating period size, \mum')
title(sprintf('Neuron %d', neuronID), 'fontsize',14)

% Adding a fit to the curve

% Trying to fit the pdf of a gamma distribution
% See curve fitting toolbox for explanation of the code
f = @(p,x) p(1)*gampdf(x/p(2),p(3),p(4));
w = 1./e;
% Infinite weights should be removed
if sum(w~=Inf)
    w(w==Inf) = max((w(w~=Inf))*5); 
else
    w(w==Inf) = 1;
end

% Specifying fitting options
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0.1 1 2 0.01],'Upper',[max(y)*3 10000 100 20]);
st_ = [1 50 10 1];
set(fo_,'Startpoint',st_);
set(fo_,'Weight',w);

ft_ = fittype('a*gampdf(x/b,c,d)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'b', 'c', 'd'});

% Doing the fitting
cf_ = fit(x,y,ft_,fo_);
p = coeffvalues(cf_);

%  Adding the fit to the plot
line(linspace(-0,max(x),400),f(p,linspace(-0,max(x),400)),...
    'color','r')

% Computing the threshold
f_thresh = @(x) p(1)*gampdf(x/p(2),p(3),p(4))-0.5;
threshold = fzero(f_thresh,mean(x),optimset('Display','off'));

end % plotSingleNeuronModulation