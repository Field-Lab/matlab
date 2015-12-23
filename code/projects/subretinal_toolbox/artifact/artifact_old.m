classdef artifact < handle
    % An artifact class that will store the average artifact shape computed
    % from a Vision .bin file. 
    % It provides with a method to remove that average trace from raw data.
    % See the accompanying processData.m file for more information.
    %
    % Version: 5.00 - 10/08/2012
    %   
    
    % Class properties
    properties
        visionPath              % path to the Vision jar archive
        rawDataFolder           % path to the folder with the raw data file
        processedFolder         % path to the folder where the average trace is saved
        startSampleExperiment   % the sample number for which experiment actually starts
        pulseDuration           % duration of the pulse
        nTrials                 % number of trials over which the artifact has to be computed
        nSkip                   % number of pulses skipped at the beginning of the data for avTrace computation
        Fs                      % sampling frequency in the original bin file - move to read-only?
        FExperiment             % Frequency of the experiment
        saveData                % Boolean, 1 if we save the artifact, 0 (default) otherwise
        useSmoothing            % 1 if you want to smooth beginning and ending of pulse, 0 otherwise
    end % public properties
    properties (SetAccess = private)
        avTrace                 % the actual artifact
        pulsePos                % the position of the pulse inside each artifact
        nElectrodes             % number of electrodes in the recording - trigger counts as 1
        nSamplesPulse           % number of samples in a pulse
        nSamplesArtifact        % number of samples in the artifact
        nSamplesSmoothing       % number of samples over which beginning and ending of pulse is smoothed
        electrodeRef            % electrode on which the maximal deflection is found
    end % read-only properties
    properties (GetAccess = private)
        dataReference           % 3-character string that indicates which file is being processed
    end % private properties
    properties (Constant)
        avTraceOffset = 80;     % Number of samples taken into account before the beginning of the pulse
        nSamplesSmoothMin = 40; % Minimum number of samples over which the data is blanked if smoothing is enabled
        trigToPulseDelay = 2;   % Delay between leading edge of the trigger and leading edge of the pulse in number of samples
    end % constant properties
    
    % Class methods
    methods
        function obj = artifact(aRawDataFolder,FExperiment,varargin)
            % Initialization of an artifact object.
            % 
            % 
            
            obj.rawDataFolder = aRawDataFolder;
            obj.FExperiment = FExperiment;

            % Setting default values for the other properties
            obj.pulseDuration = -1;
            obj.startSampleExperiment = 0;
            obj.nSkip = 10;
            obj.processedFolder = obj.rawDataFolder;
            obj.nTrials = 200;
            if isunix
                obj.visionPath = '/home/ggoetz/Research/Eclipse/110314 - Write Data V4/WriteDataFile.jar';
            else
                obj.visionPath = 'C:\Users\ggoetz\Research\Eclipse\110314 - Write Data V4\WriteDataFile.jar';
            end
            obj.saveData = 0;
            obj.useSmoothing = 0;
            experimentRun = -1;         % TODO: figure out if this is still useful
            
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
                    case 'pulseduration'
                        obj.pulseDuration = varargin{kk*2};
                    case 'ntrials'
                        obj.nTrials = varargin{kk*2};
                    case 'nskip'
                        obj.nSkip = int32(varargin{kk*2});
                    case 'processedfolder'
                        obj.processedFolder = varargin{kk*2};
                    case 'visionpath'
                        obj.visionPath = varargin{kk*2};
                    case 'startsampleexperiment'
                        obj.startSampleExperiment = int32(varargin{kk*2});
                    case 'experimentrun'
                        experimentRun = varargin{kk*2};
                    case 'savedata'
                        obj.saveData = varargin{kk*2};
                    case 'usesmoothing'
                        obj.useSmoothing = varargin{kk*2};
                    otherwise
                        err = MException('MATLAB:InvArgIn',...
                            'Unknown parameter specified');
                        throw(err);
                end
            end
            
            % Making sure vision has been linked to
            if ~exist('edu/ucsc/neurobiology/vision/io/RawDataHeader','class')
                javaaddpath(obj.visionPath);
            end
            
            % Making sure the processed folder ends by a file separator
            if obj.processedFolder(end:end)~=filesep
                obj.processedFolder = [obj.processedFolder filesep];
            end
            % Making sure it exists
            if ~exist(obj.processedFolder,'dir')
                mkdir(obj.processedFolder)
            end
            
            % If the experiment run was specified we get a datareference
            % tag
            if experimentRun >= 0
                obj.dataReference = [repmat('0',1,3-length(num2str(experimentRun))) num2str(experimentRun)];
            end
            
            % Updating the dependent properties
            obj.updateProperties();
        end % artifact
        function updateProperties(obj)
            % Should some of the public properties of the artifact be
            % changed, calling updateProperties will update the private
            % dependent properties accordingly.
            % 
            
            % Finding the missing information
            rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(obj.rawDataFolder);
            header = rawFile.getHeader();
            rawFile.close();
            
            % Sampling frequency in number of samples per ms
            obj.Fs = header.getSamplingFrequency();
            
            % Setting the number of samples in a trial
            obj.nSamplesPulse = int32(obj.Fs/obj.FExperiment);
            
            % Setting the number of samples in the artifact
            obj.nSamplesArtifact = min(max(2500,obj.pulseDuration*obj.Fs/1000*60),floor(0.8*obj.Fs/obj.FExperiment));
            
            % Setting the number of electrodes
            obj.nElectrodes = header.getNumberOfElectrodes();
            
            % Setting the smoothing width
            obj.nSamplesSmoothing = obj.nSamplesSmoothMin + ceil(obj.pulseDuration);
        end % updateProperties
        function load(obj)
            % Loads the average trace from a avArtifact.mat file saved
            % under obj.processedFolder.
            % If the .mat file is not present or if the trigger was
            % discarded in the computation of the artifact, then the
            % methods throws back an error.
            %
            % If upSampling>1 it also calls for the computation of the
            % shifted average traces (these are never saved and have to be
            % computed every time).
            % 
            
            load([obj.processedFolder 'avArtifact' obj.dataReference])
            
            % Checking that the trigger was included in the artifacts, 
            % otherwise recompute it
            try 
                lastTrace = avArtifact(:,obj.nElectrodes); %#ok<NODEF,NASGU>
            catch %#ok<CTCH>
                err = MException('MATLAB:Index out of bounds',...
                'Could not load the trigger trace: artifact has to be computed again.');
                throw(err);
            end
            
            obj.avTrace = avArtifact;
            
            % Computing the reference electrode
            obj.findElectrodeRef();
            
            % Computing the pulses positions
            obj.findPulsePos();
        end % load
        function compute(obj)
            % Computes an average artifact from the brin file in dataFolder
            % and saves that artifact in and avArtifactnnn.mat file where
            % nnn it the dataReference string.
            % 
            % If upSampling>1 it also calls for the computation of the
            % shifted average traces (these are never saved and have to be
            % computed every time).
            %
            
            DELAY_FOR_TRIG_SEARCH = 0.1;
            TARGET_TRIG_POS = 3000;
            TARGET_TRIG_POS_IN_ARTIFACT = 100;
            
            rawFile =  edu.ucsc.neurobiology.vision.io.RawDataFile(obj.rawDataFolder);
            obj.avTrace = zeros(obj.nSamplesArtifact,obj.nElectrodes,'int32');

            nn = 0;
            startSample = obj.nSkip*obj.nSamplesPulse + ...
                obj.startSampleExperiment - DELAY_FOR_TRIG_SEARCH*obj.Fs;
            % Checking that the startSample estimation is ok - times can
            % shift around
            trig_data = int32(rawFile.getData(0,startSample,obj.nSamplesPulse));
            triggerPosition = find(trig_data<(mean(trig_data)-400),1,'first');
            if triggerPosition~=TARGET_TRIG_POS
                startSample = startSample + triggerPosition - TARGET_TRIG_POS;
            end
            
            while nn<obj.nTrials
                % Getting the data
                trig_data = int32(rawFile.getData(0,startSample,obj.nSamplesPulse));

                % Checking that there was a pulse
                triggerPosition = find(trig_data(:,1)<(mean(trig_data(:,1))-400),1,'first');
                
                % Checking that the pulse has no shifted
                if triggerPosition~=TARGET_TRIG_POS
                    startSample = startSample + triggerPosition - TARGET_TRIG_POS;
                end

                % If there was one, aligning the data on the trigger pulse and
                % adding the contribution to the average trace
                if triggerPosition
                    data = int32(rawFile.getData(startSample + TARGET_TRIG_POS - TARGET_TRIG_POS_IN_ARTIFACT...
                        ,obj.nSamplesPulse));
                    if (TARGET_TRIG_POS_IN_ARTIFACT-obj.avTraceOffset+...
                            obj.nSamplesArtifact)<size(data,1)
                        obj.avTrace = obj.avTrace + data((TARGET_TRIG_POS_IN_ARTIFACT-obj.avTraceOffset+1):...
                                (TARGET_TRIG_POS_IN_ARTIFACT-obj.avTraceOffset+obj.nSamplesArtifact), :);
                    
                        % Updating the number of pulses used
                        nn = nn+1;
                    end
                end
                
                % Updating sample counter
                startSample = startSample + obj.nSamplesPulse;
            end
            
            % Normalizing the artifact
            obj.avTrace = double(obj.avTrace)/obj.nTrials;

            % Removing the dc offset from the trace
            obj.avTrace = obj.avTrace - ones(size(obj.avTrace,1),1)*obj.avTrace(1,:);
            
            % Computing the reference electrode after DC offset removal
            obj.findElectrodeRef();
            
            if obj.saveData
                avArtifact = obj.avTrace; %#ok<NASGU>
                save([obj.processedFolder 'avArtifact' obj.dataReference],'avArtifact')
            end
            
            % Computing the pulse positions
            obj.findPulsePos();
            
            % Converting to int16 type
            obj.avTrace = int16(obj.avTrace);
            
            rawFile.close();
        end % compute
        function procData = processData(obj,rawData)
            % Removes the average artifact from rawData and returns a
            % cleaned version of the data in procData.
            %
            
            [procData,posEst] = obj.removeFromData(rawData);
            if obj.useSmoothing
                procData = obj.smootheData(procData,posEst);
            end
            
        end % processData
    end % methods
    methods (Access = private)
        function findPulsePos(obj)
            % Computes the position of each artifact in the average trace.
            %
            
            % Computing the position of each pulse in the artifact
            % (by thresholding the first derivative of the pulse at 10% its 
            % maximum absolute value)
            % We look only at 1/10 of the artifact length
            avTraceGrad = abs(diff(obj.avTrace(1:(obj.avTraceOffset+floor(obj.nSamplesArtifact/10)),:,:),1,1));
            obj.pulsePos = zeros(obj.nElectrodes,1);
            for nn=1:obj.nElectrodes
                obj.pulsePos(nn) = find(avTraceGrad(:,nn)>=.10*max(avTraceGrad(:,nn)),1,'first');
            end
            obj.pulsePos = obj.pulsePos+1;
            
            % Checking that no pulse was found before the trigger 
            % (can happen for very small amplitudes)
            posTrigger = obj.pulsePos(1,:);
            
            if obj.pulsePos(obj.electrodeRef)<posTrigger(1)
                for nn=1:obj.nElectrodes
                    obj.pulsePos(nn) = posTrigger(1);
                end
            else
                for nn=1:obj.nElectrodes
                    if obj.pulsePos(nn)<posTrigger(1)
                        obj.pulsePos(nn) = obj.pulsePos(obj.electrodeRef);
                    end
                end
            end
            
        end % findPulsePos
        
        function findElectrodeRef(obj)
            % Finds the electrode with the fastest varying signal
            
            obj.electrodeRef = max(abs(diff(obj.avTrace(:,2:end,1))),[],1);
            obj.electrodeRef = find(max(obj.electrodeRef)==obj.electrodeRef,1,'first')+1;
            
        end % findElectrodeRef

        function [procData,posEst] = removeFromData(obj,rawData)
            % Actual removal method.
            % Returns the cleaned version of the data as well as posEst
            % which is an indicator of where the artifact was found inside the raw recording.. This
            % the raw recording. This value can be used in the smoothing 
            % of the data, should it be required.
            
            % Splitting artifact and trigger
            avArtifact = obj.avTrace;
            avArtifact(:,1,:) = zeros(size(obj.avTrace,1),1,size(obj.avTrace,3),'int16');

            % Trigger data
            trigger = rawData(:,1);
            
            % First check that there was a pulse
            isPulse = find(trigger<(mean(trigger)-400),1,'first');

            if isPulse
                % If there was one, then use the trigger to align and process the
                % data
                
                posEst = isPulse - obj.pulsePos(1); 

                % Sometimes the data available can be shorter than the
                % artifact, make sure it still works...
                if posEst + obj.nSamplesArtifact > size(rawData,1)
                    procData = [rawData(1:posEst,:);...
                                rawData((posEst+1):end,:) - avArtifact(1:(size(rawData,1)-posEst),:)];
                else
                    procData = [rawData(1:posEst,:);...
                                rawData((posEst+1):(posEst+obj.nSamplesArtifact),:) - avArtifact(:,:);...
                                rawData(posEst+obj.nSamplesArtifact+1:end,:)];
                end
            else
                procData = rawData;
                posEst = -1;
            end

        end % removeFromData 
        
        function procData = smootheData(obj,rawData,posEst)
            % Smoothes the beginning and ending of the pulse detected at
            % posEst by setting all the values of the recording to the 
            % average dc offset over a period of time that depends on the
            % pulse duration.
            % Smoothes only if posEst>=0, as posEst<0 means no pulse was
            % found in the data.
            %
            
            procData = rawData;
            
            if posEst>=0
                % Estimating DC offset over 100 samples of the recording
                dcOffset = mean(rawData(1:100,2:end));

                % Smoothing around artifact edges 
%                 mm=(-obj.nSamplesSmoothing):obj.nSamplesSmoothing;
                for nn=2:obj.nElectrodes
                    thisPulsePos = obj.pulsePos(nn) + posEst;
                    thisPulsePos = thisPulsePos + obj.trigToPulseDelay;
                    
                    % next line: /1000 because pulseDuration in ms and Fs in Hz
                    thisPulsePosEnd = thisPulsePos + obj.pulseDuration*obj.Fs/1000; 

                    procData(max((thisPulsePos-obj.nSamplesSmoothing),1):...
                             (thisPulsePosEnd+obj.nSamplesSmoothing),nn)...
                                   = dcOffset(1,nn-1);
%                     procData(thisPulsePos+mm,nn) = dcOffset(1,nn-1);
%                     procData(thisPulsePosEnd+mm,nn) = dcOffset(1,nn-1);
                end
            end
        end % smoothData
        
    end % private methods
    
end % classdef