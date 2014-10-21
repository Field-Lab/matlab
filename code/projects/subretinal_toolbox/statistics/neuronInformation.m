classdef neuronInformation < handle
    % Description: a class that handles neuron spikes information. 
    % Computes the PSTH and if specified the sliding PSTH of the neuron. 
    % 
    % Version: v5.01 - 04/22/2013
    %
    % Note: Spike latencies might be off by 5ms or so. To be investigated.
    %       Also, check variance calculation.
    
    % Class properties
    properties
        rawSpikeTimes           % raw spike times as returned by Vision
        ptInformation           % TTL pulses times information. Struct with time of all the pulses in the experiment.
        FExperiment             % number of light pulses per second
        stimType                % Type of stimulation: 0=pulsed IR, other=?
    end % public properties
    properties (SetAccess = private)
        neuronID                % the Vision ID for this neuron
        nExperiments            % number of different experiments analyzed 
        Fs                      % sampling frequency in the original bin file in Hertz
        spikeTimes              % spike times split between experiments
        spikeCounts             % total number of spikes in each experiment
        PSTHs                   %
        binSizes                %
        electrode               % electrode on which the neuron was found
        tags                    % sliding psth category for each experiment
        ttlsPerPulse            % Number of TTL pulses used for each PSTH
    end % read-only properties
    properties (GetAccess = private)
        nSamplesTrial           % number of samples per pulse
    end % private properties
    properties (Constant)
        defaultFExperiment = 2  % for spontaneous/CW vis data, the PSTH has to be computed,
                                % and it is done using this dummy frequency
        defaultBinSize = 100    % default bin size for the PSTH
    end % constant properties
    
    % Class methods
    methods
        function obj = neuronInformation(varargin)
            if nargin==0
                obj.neuronID = -1;
            else
                % Setting default values
                obj.Fs = 20000;
                obj.ttlsPerPulse = 1;
                
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
                            ['Unexpected property specified in position ' num2str(kk)]);
                        throw(err);
                    end

                    switch lower(varargin{kk*2-1})
                        case 'neuronid'
                            obj.neuronID = varargin{kk*2};
                        case 'rawspiketimes'
                            obj.rawSpikeTimes = double(varargin{kk*2});
                        case 'pulsetimeinformation'
                            obj.ptInformation = varargin{kk*2};
                        case 'fexperiment'
                            obj.FExperiment = varargin{kk*2};
                        case 'fs'
                            obj.Fs = varargin{kk*2};
                        case 'electrode'
                            obj.electrode = varargin{kk*2};
                        case 'stimulationtype'
                            obj.stimType = varargin{kk*2};
                        case 'ttlsperpulse'
                            obj.ttlsPerPulse = varargin{kk*2};
                        otherwise
                            err = MException('MATLAB:InvArgIn',...
                                'Unknown parameter specified');
                            throw(err);
                    end
                end

                obj.nExperiments = numel(obj.ptInformation);
                
                % If no information about the stimulation type, we assume
                % pulsed IR
                if isempty(obj.stimType)
                    obj.stimType = zeros(obj.nExperiments,1);
                end
                
                % Then filling in FExperiment if it was not specified
                if isempty(obj.FExperiment)
                    obj.FExperiment = obj.defaultFExperiment*ones(obj.nExperiments,1);
                else if length(obj.FExperiment)==1
                        obj.FExperiment = obj.FExperiment*ones(obj.nExperiments,1);
                    end
                end
                % For the non pulsed IR stimulation, we replace the
                % meaningless FExperiment (which can be negative) by the
                % default value, it is useful for computing PSTH and will
                % make things go wrong if not done.
                for kk=1:obj.nExperiments
                    if obj.stimType(kk)~=0
                        obj.FExperiment(kk) = obj.defaultFExperiment;
                    end
                end
                    
                obj.nSamplesTrial = round(obj.Fs*(1./obj.FExperiment)*obj.ttlsPerPulse);
                
                obj.splitRawSpikeTimes();
            end
        end % neuronData
        
        function computeFiringInformation(obj,binSize)
            % Computes the PSTH, sliding PSTH and categorizes each
            % experiment depending on the look of the sliding PSTH
            % If you want to override the optimal bin size detection
            % happening in the computation, it is possible to provide 
            % the function with a bin size that will be used instead.
            % This binsize should be specified in number of samples and 
            % not in milliseconds.
            % 
            
            if nargin==0
                overrideBinSize = false;
            else
                overrideBinSize = true;
            end
            
            obj.PSTHs = struct('experimentID',{},...
                               'data',{},...
                               'variances',{});
            obj.binSizes = struct('experimentID',{},...
                                  'data',{});
            obj.tags = struct('experimentID',{},...
                              'data',{});
            
            for kk=1:obj.nExperiments
                obj.PSTHs(kk).experimentID = kk;
                obj.binSizes(kk).experimentID = kk;
                if overrideBinSize
                    obj.computePSTH(kk,binSize);
                else
                    obj.computePSTH(kk);
                end
            end
            
        end % computeFiringInformation
        
        function saveFigures(obj,dataFolder, varargin)
            % Plots the sliding PSTH computed and saves them in dataFolder.
            % it is possible to specify the file extension (by default,
            % .eps) and change it to any file format supported by Matlab.
            %
            % Optional parameters: imageFormat, displayImage, saveImage,
            % combineExperiments, nameTag
            %
            % combineExperiments should be a structure with the following
            % fields:
            %   - varyingParam: string with the name of the parameter
            %   that's varying over all the experiments we wish to plot
            %   (for example 'pulse width', 'power')
            %   - unitParam: unit of the varying param
            %   - expID: the ID of the experiments we're combining
            %
            % nameTag is an optional tag that will go to the end of the
            % file name when it is saved (useful for example when several
            % combineExperiment parameters are called successively, without
            % a change in varying parameter)
            %
            
            if dataFolder(end:end)~=filesep
                dataFolder = [dataFolder filesep];
            end
            imageFormat = 'eps';
            
            % Setting default values for optional input parameters
            displayImage = 'off';
            combineExperiments = struct('varyingParam',{},'unitParam',{},'valuesParam',{},'expID',{});
            saveImage = true;
            nameTag = '';
            % Reading the optional input arguments
            nbin = length(varargin);
            for kk=1:(nbin/2)
                if ~ischar(varargin{kk*2-1})
                    err = MException('MATLAB:InvArgIn',...
                        ['Unexpected property specified in position ' num2str(kk)]);
                    throw(err);
                end

                switch lower(varargin{kk*2-1})
                    case 'imageformat'
                        imageFormat = varargin{kk*2};
                    case 'displayimage'
                        displayImage = varargin{kk*2};
                    case 'saveimage'
                        saveImage = varargin{kk*2};
                    case 'nametag'
                        nameTag = ['_' varargin{kk*2}];
                    case 'combineexperiments'
                        combineExperiments = varargin{kk*2};
                    otherwise
                        err = MException('MATLAB:InvArgIn',...
                            'Unknown parameter specified');
                        throw(err);
                end
            end
            
            if ~exist([dataFolder num2str(obj.neuronID)],'dir')
                mkdir([dataFolder num2str(obj.neuronID)]);
            end
            
            fh = figure; clf; set(fh,'color','white','visible',displayImage);
            
            if isempty(combineExperiments)
                for kk=1:obj.nExperiments
                    % Plotting the data
                    stairs((0:obj.binSizes(kk).data:obj.nSamplesTrial(kk))/...
                        obj.Fs,obj.PSTHs(kk).data)
                    
                    % Axis and labels
%                     axis([0 obj.nSamplesTrial(kk)/obj.Fs 0 100])
                    title(['Neuron ID: ' num2str(obj.neuronID)...
                           '; Experiment run: ' num2str(kk)...
                            '; Electrode: ' num2str(obj.electrode)...
                            '; Number of spikes: ' num2str(obj.spikeCounts(kk).data)])
                    xlabel('Time after stimulation pulse, ms')
                    ylabel('#spikes/trial/bin');
                    hold off

                    if saveImage
                        saveas(fh,[dataFolder num2str(obj.neuronID) filesep...
                            'neuron' num2str(obj.neuronID) '_run' num2str(kk)...
                            nameTag'],imageFormat)
                    end

                    clf;
                end
            else
                % Combining the data
                numSpikes = 0;
                xAxis = [];
                plotData = [];
                
                % In case the frequency was changed between experiments, we
                % first check how long the longest PSTH will be
                maxLength = 0;
                for kk=1:length(combineExperiments.expID)
                    thisXDataStairs = (0:obj.binSizes(combineExperiments.expID(kk)).data...
                        :obj.nSamplesTrial(combineExperiments.expID(kk)))/obj.Fs;
                    if length(thisXDataStairs)>maxLength
                        xDataStairs = thisXDataStairs;
                        maxLength = length(thisXDataStairs);
                    end
                end
                
                for kk=1:length(combineExperiments.expID)
                    % Getting the points we'll later plot.
                    % The x axis data is the same for all experiments so
                    % that different frequencies can be combined
                    % If there was not enough y data (should happen when
                    % combining frequencies) we append 0s at the end of y
                    yDataStairs = [obj.PSTHs(combineExperiments.expID(kk)).data...
                        zeros(1, maxLength-length(...
                        obj.PSTHs(combineExperiments.expID(kk)).data))];

                    % Getting the stairs plot
                    [xAxis, plotData(kk,:)] = stairs(xDataStairs,yDataStairs);
                    
                    numSpikes = numSpikes + obj.spikeCounts(combineExperiments.expID(kk)).data;
                end
                
                % Plotting the data
                plot(xAxis,plotData)
                
                % Axis
                axis([0 max(xAxis) 0 max(max(plotData))*1.1+10^(-6)])
                % Title, axis label
                title(['Neuron ID: ' num2str(obj.neuronID)...
                       '; Varying '  combineExperiments.varyingParam...
                       '; Electrode: ' num2str(obj.electrode)...
                       '; Total number of spikes: ' num2str(numSpikes)])
                xlabel('Time after stimulation pulse, s')
                ylabel('Number of spikes per trial and per bin');
                % Plot legend: creating the string matrix
                legendStr = [num2str(combineExperiments.valuesParam(1),'%1.1f')...
                    ' ' combineExperiments.unitParam];
                
                for kk=2:length(combineExperiments.expID)
                    currentStr = [num2str(combineExperiments.valuesParam(kk),'%1.1f')...
                        ' ' combineExperiments.unitParam];
                    
                    if size(legendStr,2)<size(currentStr,2)
                        fillStr = repmat(' ',size(legendStr,1),...
                            size(currentStr,2)-size(legendStr,2));
                        
                        legendStr = [fillStr legendStr];
                    else
                        if size(legendStr,2)>size(currentStr,2)
                            fillStr = repmat(' ',1,size(legendStr,2)-size(currentStr,2));
                            currentStr = [fillStr currentStr];
                        end
                    end
                    
                    legendStr = [legendStr; currentStr];
                end
                % Finally adding the legend
                legend(legendStr);
                
                if strcmp(displayImage,'on')
                    pause
                end
                
                if saveImage
                    saveas(fh,[dataFolder num2str(obj.neuronID) filesep...
                        'neuron' num2str(obj.neuronID) '_' ...
                        combineExperiments.varyingParam '_combined' ...
                        nameTag],imageFormat)
                end
                
                clf;
            end
            close(fh)
        end % createFigures
        
        function save(obj,dataFolder)
            % Saves the neuronInformation instance in the specified data
            % folder
            
            if dataFolder(end:end)~=filesep
                dataFolder = [dataFolder filesep];
            end
            save([dataFolder 'neuron' num2str(obj.neuronID) '.mat'],'obj');
        end % save
        
        function allPSTHs = getPSTHs(obj)
            % Returns a matrix with all the PSTHs data for that particular
            % neuron.
            
            nPoints = length(obj.PSTHs(1).data);
            allPSTHs = zeros(obj.nExperiments,nPoints);
            for kk=1:obj.nExperiments
                allPSTHs(kk,:) = obj.PSTHs(kk).data;
            end
        end
        
        function theElectrode = getElectrode(obj)
            % Returns the electrode number on which the neuron was found
            theElectrode = obj.electrode;
        end % getElectrode
        
        function nSpikes = getSpikeCounts(obj)
            % Returns the number of spikes for this neuron for each
            % experiment
            
            nSpikes = zeros(obj.nExperiments,1);
            for kk=1:obj.nExperiments
                nSpikes(kk) = obj.spikeCounts(kk).data;
            end
        end % getNumberOfSpikes
        
        function theID = getNeuronID(obj)
            % Returns the neuron ID of that neuron
            theID = obj.neuronID;
        end % getNeuronID
        
        function nExpRun = getNumberOfExperiments(obj)
            % Returns the number of experiments that were analyzed for this
            % neuron
            nExpRun = obj.nExperiments;
        end % getNumberOfExperiments
        
    end % public methods
    
    methods (Access = private)
        
        function splitRawSpikeTimes(obj)
            % Computes the spike times for each experiment
            % t=0 corresponds to the time of the first pulse of light in
            % each experiment, so t<0 is possible and should not be
            % alarming.
            
            obj.spikeTimes = struct('experimentID',{},...
                                    'data',{});
            obj.spikeCounts = struct('experimentID',{},...
                                    'data',{});
                                
            for kk=1:obj.nExperiments
                obj.spikeTimes(kk).experimentID = kk;
                obj.spikeCounts(kk).experimentID = kk;

                startSample = obj.ptInformation(kk).startSample;
                endSample = startSample + obj.ptInformation(kk).nSamplesExperiment;
                
                % Keep only the spikes that happened between startSample
                % and endSample.
                obj.spikeTimes(kk).data = obj.rawSpikeTimes(logical(...
                                            (obj.rawSpikeTimes>startSample)...
                                            .*(obj.rawSpikeTimes<=endSample)));
                obj.spikeTimes(kk).data = obj.spikeTimes(kk).data - startSample;
                
                obj.spikeCounts(kk).data = numel(obj.spikeTimes(kk).data);
            end
        end % splitRawSpikeTimes
        
        function wrappedSpikeTimes =  wrapSpikeTimes(obj,spikeTimes,experimentID)
            % For a given experiment, returns a matrix with the spike times
            % that happened during that experiment. 
            % 
            % If spikeTimes is a matrix with spike times, the result will
            % be a matrix of time elapsed (in samples) between the last
            % stimulation pulse and the spike times. 
            % 
            % The resulting wrappedSpikeTimes matrix is sorted and no value
            % in it will be greater than the number of samples in a given
            % trial for experiment experimentID.
            
            currentPulseTimes = obj.ptInformation(experimentID).TTLTimes;
            currentPulseTimes = currentPulseTimes(1:obj.ttlsPerPulse:end);

            if length(currentPulseTimes)==1
                % If there is no pulse time in between: it's not
                % stimulation data and the wrapping can be basic
                nSamplesPulse = obj.nSamplesTrial(experimentID);
                wrappedSpikeTimes = sort(mod(spikeTimes,nSamplesPulse));
            else
                % Otherwise, have to take into account each pulse time
                wrappedSpikeTimes = zeros(size(spikeTimes));
                for kk=1:length(spikeTimes)
                    latestPulse = find(currentPulseTimes<spikeTimes(kk),1,'last');
                    latestPulse = currentPulseTimes(latestPulse);
                    wrappedSpikeTimes(kk) = spikeTimes(kk) - latestPulse;
                end
            end
            
        end % wrapSpikeTimes
        
        function computePSTH(obj,experimentID,binSize)
            % Computes the optimal bin size for the spike times
            % distribution in the experiment analyzed, and creates the
            % corresponding PSTH, if not provided with a bin size.
            % If provided with a bin size, computes the corresponding PSTH.
            % 
            
            currentPulseTimes = obj.ptInformation(experimentID).TTLTimes;
            currentPulseTimes = currentPulseTimes(1:obj.ttlsPerPulse:end);
            nPulses = length(currentPulseTimes);

            nSamplesPulse = obj.nSamplesTrial(experimentID);
            
            % The way the PSTH is built depends on the stimulus type 
            switch obj.stimType(experimentID)
                % Pulsed IR or pulsed visible case
                case {0,1}
                    if nargin==1
                        currentSpikeTimes = obj.wrapSpikeTimes(...
                            obj.spikeTimes(experimentID).data,experimentID);

                        % We don't compute the variance
                        obj.PSTHs(experimentID).variances = NaN;

                        % Computing the PSTH
                        C = Inf;    % cost function
                        for binSize = floor(logspace(log10(25),log10(500),100))
                            % Binning the spike data
                            PSTH = histc(currentSpikeTimes,0:binSize:nSamplesPulse);
                            % Computing the corresponding cost function
                            kbar = mean(PSTH);
                            nu = var(PSTH);
                            CbinSize = (2*kbar - nu)/(nPulses*binSize/obj.Fs)^2;
                            % If the estimate is better then we change the PTSH and binsize
                            if CbinSize < C
                                C = CbinSize;
                                % Storing the scaled version of the PTSH
                                thePSTH = PSTH/nPulses; % to get firing frequency replace nPulses by (binSize/20000*nPulses);
                                thebinSize = binSize;
                            end
                        end
                    else % we also get the variance
                        cSubPSTH = zeros(nPulses,length(0:binSize:nSamplesPulse));
                        experimentData = obj.spikeTimes(experimentID).data;
                        for ll=1:nPulses
                            % Getting the spike times for that one pulse
                            if ll<nPulses
                                startSample = double(currentPulseTimes(ll));
                                endSample = double(currentPulseTimes(ll+1));
                            else
                                startSample = double(currentPulseTimes(ll));
                                endSample = startSample + nSamplesPulse;
                            end

                            currentData = experimentData(logical((experimentData...
                                >=startSample).*(experimentData<=endSample)));

                            % Putting them in a sub-histogram;
                            if ~isempty(currentData)
                                % No need for the wrapSpikesTimes method in this case:
                                % we just have to remove startSample from the data
                                currentData = currentData - startSample;
                                
                                cSubPSTH(ll,:) = histc(currentData,0:binSize:nSamplesPulse);
                            end
                        end

                        % Binning the spike data
                        thePSTH = sum(cSubPSTH,1)/nPulses; % to get firing frequency replace nPulses by (binSize/20000*nPulses)
                        thebinSize = binSize;

                        obj.PSTHs(experimentID).variances = var(cSubPSTH,0,1);
                        obj.PSTHs(experimentID).numberOfPulses = nPulses;
                        obj.PSTHs(experimentID).subPSTH = cSubPSTH;
                    end
                
                % CW visible case/spontaneous
                case 2
                    if nargin==1
                        binSize = obj.defaultBinSize;
                    end
                    currentSpikeTimes = obj.wrapSpikeTimes(...
                            obj.spikeTimes(experimentID).data,experimentID);
                    
                    % No variance
                    obj.PSTHs(experimentID).variance = NaN;
                    
                    % PSTH + scaling factor and binSize
                    thePSTH = histc(currentSpikeTimes,0:binSize:nSamplesPulse);
                    thePSTH = thePSTH/(obj.nSamplesExperiment(experimentID)/nSamplesPulse);
                    thebinSize = binSize;
            end
            
            obj.PSTHs(experimentID).data = thePSTH;
            obj.binSizes(experimentID).data = thebinSize;
        end % computePSTH
        
    end % private methods
    
end % classdef

