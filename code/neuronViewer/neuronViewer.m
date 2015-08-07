classdef neuronViewer < handle
    %NEURONVIEWER Simple class for cluster visualization
    %   Constructor loads a .neurons-raw and a .prj file
    %   If available, it also loads a .model file
    %
    % Author -- Vincent Deo -- Stanford University -- August 7, 2015
    
    
    properties (SetAccess = immutable, GetAccess = public)
        prjFilePath
        modelFilePath
        neuronFilePath
        
        prjFile
        modelFile
        neuronFile
        
        nElectrodes
    end
    
    properties (Access = private)
        projSpikes
        
        numNeurons
        neuronIDs
        neuronEls
        neuronSpikes
        
        neurSubset
        
        elLoaded
        spikesEl
        prjEl
        clusters
        modelMeans
        modelCovariances
    end
    
    methods
        % Constructor
        %
        % Input: path to analysis folder
        function obj = neuronViewer(analysisPath,varargin)
            
            % Argument check
            validateattributes(analysisPath,{'char'},{},'','Analysis Path',1);
            if ~ (exist(analysisPath,'file') == 7)
                throw(MException('','neuronViewer:neuronViewer - Analysis folder not found'));
            end
            
            if nargin > 1
                narginchk(4,4);
                validateattributes(varargin{1},{'char'},{},'','projections file',2);
                validateattributes(varargin{2},{'char'},{},'','model file',3);
                validateattributes(varargin{3},{'char'},{},'','neurons-raw file',4);
            end
            
            % Add custom classes and vision to the path
            javaaddpath('/Volumes/Lab/Users/vdeo/matlab/code/neuronViewer/');
            if ~(exist('edu.ucsc.neurobiology.vision.io.NeuronFile','class'));
                javaaddpath('/Volumes/Lab/Development/vision7/Vision.jar','-end');
                % javaaddpath('C:\Users\Vincent\Documents\EJGroup\mvision\vision\Vision.jar','-end');
            end
            import edu.ucsc.neurobiology.vision.io.*
            
            % Look for a .model, a .prj and a .neurons raw
            if nargin == 1 || (nargin == 4 && numel(varargin{1}) == 0)
                files = dir([analysisPath,filesep,'*.prj']);
            else
                files = dir([analysisPath,filesep,varargin{1}]);
            end
            if numel(files) == 1
                obj.prjFilePath = [analysisPath,filesep,files(1).name];
                obj.prjFile = edu.ucsc.neurobiology.vision.io.ProjectionsFile(obj.prjFilePath);
            else if numel(files) == 0
                    throw(MException('','neuronViewer:neuronViewer - No projections file found'));
                else
                    throw(MException('','neuronViewer:neuronViewer - Multiple projections files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,modelFile,neuronsRawFile)'));
                end
            end
            
            if nargin == 1 || (nargin == 4 && numel(varargin{3}) == 0)
                files = dir([analysisPath,filesep,'*.neurons-raw']);
            else
                files = dir([analysisPath,filesep,varargin{3}]);
            end
            if numel(files) == 1
                obj.neuronFilePath = [analysisPath,filesep,files(1).name];
                obj.neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(obj.neuronFilePath);
            else if numel(files) == 0
                    throw(MException('','neuronViewer:neuronViewer - No neurons-raw file found'));
                else
                    throw(MException('','neuronViewer:neuronViewer - Multiple neurons-raw files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,modelFile,neuronsRawFile)'));
                end
            end
            
            if nargin == 1 || (nargin == 4 && numel(varargin{2}) == 0)
                files = dir([analysisPath,filesep,'*.model']);
            else
                files = dir([analysisPath,filesep,varargin{2}]);
            end
            if numel(files) == 1
                obj.modelFilePath = [analysisPath,filesep,files(1).name];
                obj.modelFile = edu.ucsc.neurobiology.vision.io.ClusteringModelFile(obj.modelFilePath);
            else if numel(files) == 0
                    disp('neuronViewer:neuronViewer - No model file found, skipping');
                    obj.modelFilePath = [];
                    obj.modelFile = [];
                else
                    throw(MException('','neuronViewer:neuronViewer - Multiple model files found\nPlease use extended constructor call\nneuronViewer(analysisPath,projectionsFile,modelFile,neuronsRawFile)'));
                end
            end
            
            obj.numNeurons = obj.neuronFile.getNumberOfNeurons();
            obj.neuronIDs = obj.neuronFile.getIDList();
            obj.neuronEls = zeros(numel(obj.numNeurons),1);
            
            for neurIndex = 1:obj.numNeurons
                obj.neuronEls(neurIndex) = obj.neuronFile.getElectrode(obj.neuronIDs(neurIndex));
            end
            
            obj.elLoaded = -1; % JAVA numbering
            
            h = obj.prjFile.getHeader();
            id = h.arrayID();
            elMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(id);
            obj.nElectrodes = elMap.getNumberOfElectrodes;
            
        end % Constructor
        
        % Destructor
        function delete(obj)
            import edu.ucsc.neurobiology.vision.io.*
            obj.prjFile.close();
            if numel(obj.modelFile) > 0
                obj.modelFile.close();
            end
            obj.neuronFile.close();
        end % Destructor
        
        % Draw clusters
        % Syntax
        % obj.drawClusters(electrode)
        % obj.drawClusters(electrode, number plot points)
        % obj.drawClusters(electrode, number plot points, plot boxes y/n)
        % obj.drawClusters(electrode, 0, plot boxes) for no subsampling
        function drawClusters(obj, el, varargin)
            narginchk(2,4)
            validateattributes(el,{'numeric'},{'scalar','integer','>',0,'<',obj.nElectrodes},'','electrode index',1);
            if nargin >= 3
                validateattributes(varargin{1},{'numeric'},{'scalar','integer','>=',0},'','subset size',2);
            end
            if nargin >= 4
                validateattributes(varargin{2},{'numeric','logical'},{'scalar','integer','>=',0,'<=',1},'','plot boxes',3);
            end
            
            import edu.ucsc.neurobiology.vision.io.*
            
            if el ~= obj.elLoaded
                % Load neurons
                obj.neurSubset = find(obj.neuronEls == el);
                obj.neuronSpikes = cell(numel(obj.neurSubset),1);
                for i = 1:numel(obj.neurSubset)
                    obj.neuronSpikes{i} = obj.neuronFile.getSpikeTimes(obj.neuronIDs(obj.neurSubset(i)));
                end
                
                % Load projections
                prjWrap = obj.prjFile.readProjections(el);
                cut = find(prjWrap.times,1,'last');
                obj.spikesEl = prjWrap.times(1:cut);
                obj.prjEl = prjWrap.data(:,1:cut)';
                
                obj.clusters = zeros(cut,1);
                
                % Make cluster indexes
                for i = 1:numel(obj.neurSubset)
                    [~,ind,~] = intersect(obj.spikesEl,obj.neuronSpikes{i});
                    obj.clusters(ind) = i;
                end
                
                % Load model boxes
                if numel(obj.modelFile) > 0
                    model = obj.modelFile.getNeuronExtraction(el);
                    if numel(model) > 0
                        obj.modelMeans = model.means;
                        obj.modelCovariances = model.covariances;
                    else
                        obj.modelMeans = [];
                        obj.modelCovariances = [];
                    end
                end
            end
            
            obj.elLoaded = el;
            
            % subsampling
            if nargin >= 3 && varargin{1} < size(obj.prjEl,1) && varargin{1} > 0
                subsample = randsample(size(obj.prjEl,1), varargin{1});
                scatter3(obj.prjEl(subsample,1),obj.prjEl(subsample,2),obj.prjEl(subsample,3),9,obj.clusters(subsample));
            else
                scatter3(obj.prjEl(:,1),obj.prjEl(:,2),obj.prjEl(:,3),9,obj.clusters);
            end
            
            IDs = num2cell(obj.neuronIDs(obj.neurSubset));
            caxis([0 numel(obj.neurSubset)]);
            colorbar('ticklabels',[{'Outliers'};IDs]);
            colormap hsv
            
            % Drawing model boxes
            if (numel(obj.modelFile) > 0 && nargin < 4) || (numel(obj.modelFile) > 0 && nargin == 4 && varargin{2})
                hold on
                kSig = 1;
                lineMat = [1,-1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,-1,1;...
                    -1,-1,1,1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1;...
                    -1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,1,-1,1,-1];
                
                C = get(gca,'CLim');
                CM = colormap;
                for i = 1:size(obj.modelMeans,1)                    
                    Ci = floor((i-C(1))/(C(2)-C(1))*size(CM,1));
                    cube = bsxfun(@plus,obj.modelMeans(i,1:3)',kSig * bsxfun(@times,sqrt(obj.modelCovariances(i,1:3)'),lineMat));
                    % plot3(cube(1,:),cube(2,:),cube(3,:),'color',CM(Ci,:),'linewidth',2);
                    plot3(cube(1,:),cube(2,:),cube(3,:),'k-','linewidth',2);
                end
                axis tight
                hold off
            end
            
        end %drawClusters
        
        function flush(obj)
            % Do nothing so far
        end % flush
        
    end % Methods
end

