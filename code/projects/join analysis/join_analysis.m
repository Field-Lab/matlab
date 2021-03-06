% identify which spike sorting parameter set found which neurons
%
%
% 11-08 tamachado and jgauthier
%

% choose dataset
switch 6
    case 6
        main = '/Volumes/War/Analysis/Machado/2005-04-26-0/data002/params';
        fo = {'4.0-5-15-1-5-30-3-3-8',... %this is the join order
            '4.0-5-15-1-5-30-5-3-8',...
            '4.0-5-15-1-5-30-7-3-8',...
            '4.0-5-15-1-5-40-3-3-8',...
            '4.0-5-15-1-5-40-5-3-8',...
            '4.0-5-15-1-5-40-7-3-8',...
            '4.0-5-15-1-5-50-3-3-8',...
            '4.0-5-15-1-5-50-5-3-8',...
            '4.0-5-15-1-5-50-7-3-8',...
            '4.0-5-15-1-5-30-3-3-15',...
            '4.0-5-15-1-5-30-5-3-15',...
            '4.0-5-15-1-5-30-7-3-15',...
            '4.0-5-15-1-5-40-3-3-15',...
            '4.0-5-15-1-5-40-5-3-15',...
            '4.0-5-15-1-5-40-7-3-15',...
            '4.0-5-15-1-5-50-3-3-15',...
            '4.0-5-15-1-5-50-5-3-15',...
            '4.0-5-15-1-5-50-7-3-15'};
        classFilePath = '/Volumes/Blues/Data/Clare/Tim/2005-04-26-0/data002/joined/joined.classification.txt';
        CELL = 1; %fo is a cell array
    case 5
        main = '/Volumes/Bore/Analysis.noindex/Machado/2008-11-10-0/data005/params';
        fo = {'3.0-5-15-1-5-30-3-3-8',... %this is the join order
            '3.0-5-15-1-5-30-5-3-8',...
            '3.0-5-15-1-5-30-7-3-8',...
            '3.0-5-15-1-5-40-3-3-8',...
            '3.0-5-15-1-5-40-5-3-8',...
            '3.0-5-15-1-5-40-7-3-8',...
            '3.0-5-15-1-5-50-3-3-8',...
            '3.0-5-15-1-5-50-5-3-8',...
            '3.0-5-15-1-5-50-7-3-8',...
            '3.0-10-12-1-5-40-3-2-8',...
            '4.0-5-15-1-5-30-3-3-8',...
            '4.0-5-15-1-5-30-5-3-8',...
            '4.0-5-15-1-5-30-7-3-8',...
            '4.0-5-15-1-5-40-3-3-8',...
            '4.0-5-15-1-5-40-5-3-8',...
            '4.0-5-15-1-5-40-7-3-8',...
            '4.0-5-15-1-5-50-3-3-8',...
            '4.0-5-15-1-5-50-5-3-8',...
            '4.0-5-15-1-5-50-7-3-8'};
        classFilePath = '/Volumes/Cumberland/Analysis.noindex/Clare/Tim/2008-11-10-0/joined/joined.classification.txt';
        CELL = 1; %fo is a cell array
    case 4
        main = '/Volumes/Terrapin/Analysis/Machado/2005-04-26-0/data009/params';
        fo = ['4.0-10-12-1-5-30-5-2-8'
            '4.0-10-12-1-5-30-7-2-8'
            '4.0-10-12-1-5-40-3-2-8'
            '4.0-10-12-1-5-40-5-2-8'
            '4.0-10-12-1-5-40-7-2-8'
            '4.0-10-12-1-5-50-3-2-8'
            '4.0-10-12-1-5-50-5-2-8'
            '4.0-10-12-1-5-50-7-2-8'
            '5.0-10-12-1-5-30-3-2-8'
            '5.0-10-12-1-5-30-5-2-8'
            '5.0-10-12-1-5-30-7-2-8'
            '5.0-10-12-1-5-40-3-2-8'
            '5.0-10-12-1-5-40-5-2-8'
            '5.0-10-12-1-5-40-7-2-8'
            '5.0-10-12-1-5-50-3-2-8'
            '5.0-10-12-1-5-50-5-2-8'
            '5.0-10-12-1-5-50-7-2-8'];
        classFilePath = '/Volumes/Cumberland/Analysis.noindex/Clare/Tim/2005-04-26-0/joined/joined/joined-classification.txt';
        CELL = 0; %fo is a not a cell it is an array
    case 1
        main = '/Volumes/Terrapin/Analysis/Machado/2008-04-30-2/data005/params';
        fo = ['4.0-10-12-1-5-30-3-2-8'
            '4.0-10-12-1-5-50-7-5-8'
            '5.0-10-12-1-5-40-5-7-8'
            '4.0-10-12-1-5-30-3-5-8'
            '4.0-10-12-1-5-50-7-7-8'
            '5.0-10-12-1-5-50-7-3-8'
            '4.0-10-12-1-5-30-3-7-8'
            '5.0-10-12-1-5-30-3-3-8'
            '5.0-10-12-1-5-50-7-5-8'
            '4.0-10-12-1-5-40-5-3-8'
            '5.0-10-12-1-5-30-3-5-8'
            '5.0-10-12-1-5-50-7-7-8'
            '4.0-10-12-1-5-40-5-5-8'
            '5.0-10-12-1-5-30-3-7-8'
            '4.0-10-12-1-5-40-5-7-8'
            '5.0-10-12-1-5-40-5-3-8'
            '4.0-10-12-1-5-50-7-3-8'
            '5.0-10-12-1-5-40-5-5-8'];
        classFilePath = '/Volumes/Blues/Data/Clare/Tim/512-params-union/512-params-union.classification.txt';
        CELL = 0;
    case 2
        main = '/Volumes/Sunshine/Analysis/Machado/2008-08-27-5/data004/params';
        fo = ['4.0-10-12-1-5-30-3-2-8'
            '4.0-10-12-1-5-50-5-2-8'
            '5.0-10-12-1-5-40-7-2-8'
            '4.0-10-12-1-5-30-5-2-8'
            '4.0-10-12-1-5-50-7-2-8'
            '5.0-10-12-1-5-50-3-2-8'
            '4.0-10-12-1-5-30-7-2-8'
            '5.0-10-12-1-5-30-3-2-8'
            '5.0-10-12-1-5-50-5-2-8'
            '4.0-10-12-1-5-40-3-2-8'
            '5.0-10-12-1-5-30-5-2-8'
            '5.0-10-12-1-5-50-7-2-8'
            '4.0-10-12-1-5-40-5-2-8'
            '5.0-10-12-1-5-30-7-2-8'
            '4.0-10-12-1-5-40-7-2-8'
            '5.0-10-12-1-5-40-3-2-8'
            '4.0-10-12-1-5-50-3-2-8'
            '5.0-10-12-1-5-40-5-2-8'];
        classFilePath = '/Volumes/Blues/Analysis/Machado/519-params-union/joined.classification.txt';
        CELL = 0;
    case 3
        main = '/Volumes/Alligator/Analysis/Gauthier/2008-04-30-1/params';
        fo = {'3.0-5-15-0-5-30-3-2-8',...
            '4.0-5-15-0-5-30-3-2-8',...
            '5.0-5-15-0-5-30-3-2-8',...
            '3.0-10-15-0-5-30-3-2-8',...
            '4.0-10-15-0-5-30-3-2-8',...
            '5.0-10-15-0-5-30-3-2-8',...
            '3.0-5-15-1-5-30-3-2-8',...
            '4.0-5-15-1-5-30-3-2-8',...
            '5.0-5-15-1-5-30-3-2-8',...
            '3.0-10-15-1-5-30-3-2-8',...
            '4.0-10-15-1-5-30-3-2-8',...
            '5.0-10-15-1-5-30-3-2-8',...
            '3.0-5-15-2-5-30-3-2-8',...
            '4.0-5-15-2-5-30-3-2-8',...
            '5.0-5-15-2-5-30-3-2-8',...
            '3.0-10-15-2-5-30-3-2-8',...
            '4.0-10-15-2-5-30-3-2-8',...
            '5.0-10-15-2-5-30-3-2-8',...
            '3.0-5-15-0-5-40-3-2-8',...
            '4.0-5-15-0-5-40-3-2-8',...
            '5.0-5-15-0-5-40-3-2-8',...
            '3.0-10-15-0-5-40-3-2-8',...
            '4.0-10-15-0-5-40-3-2-8',...
            '5.0-10-15-0-5-40-3-2-8',...
            '3.0-5-15-1-5-40-3-2-8',...
            '4.0-5-15-1-5-40-3-2-8',...
            '5.0-5-15-1-5-40-3-2-8',...
            '3.0-10-15-1-5-40-3-2-8',...
            '4.0-10-15-1-5-40-3-2-8',...
            '5.0-10-15-1-5-40-3-2-8',...
            '3.0-5-15-2-5-40-3-2-8',...
            '4.0-5-15-2-5-40-3-2-8',...
            '5.0-5-15-2-5-40-3-2-8',...
            '3.0-10-15-2-5-40-3-2-8',...
            '4.0-10-15-2-5-40-3-2-8',...
            '5.0-10-15-2-5-40-3-2-8',...
            '3.0-5-15-0-5-30-5-2-8',...
            '4.0-5-15-0-5-30-5-2-8',...
            '5.0-5-15-0-5-30-5-2-8',...
            '3.0-10-15-0-5-30-5-2-8',...
            '4.0-10-15-0-5-30-5-2-8',...
            '5.0-10-15-0-5-30-5-2-8',...
            '3.0-5-15-1-5-30-5-2-8',...
            '4.0-5-15-1-5-30-5-2-8',...
            '5.0-5-15-1-5-30-5-2-8',...
            '3.0-10-15-1-5-30-5-2-8',...
            '4.0-10-15-1-5-30-5-2-8',...
            '5.0-10-15-1-5-30-5-2-8',...
            '3.0-5-15-2-5-30-5-2-8',...
            '4.0-5-15-2-5-30-5-2-8',...
            '5.0-5-15-2-5-30-5-2-8',...
            '3.0-10-15-2-5-30-5-2-8',...
            '4.0-10-15-2-5-30-5-2-8',...
            '5.0-10-15-2-5-30-5-2-8',...
            '3.0-5-15-0-5-40-5-2-8',...
            '4.0-5-15-0-5-40-5-2-8',...
            '5.0-5-15-0-5-40-5-2-8',...
            '3.0-10-15-0-5-40-5-2-8',...
            '4.0-10-15-0-5-40-5-2-8',...
            '5.0-10-15-0-5-40-5-2-8',...
            '3.0-5-15-1-5-40-5-2-8',...
            '4.0-5-15-1-5-40-5-2-8',...
            '5.0-5-15-1-5-40-5-2-8',...
            '3.0-10-15-1-5-40-5-2-8',...
            '4.0-10-15-1-5-40-5-2-8',...
            '5.0-10-15-1-5-40-5-2-8',...
            '3.0-5-15-2-5-40-5-2-8',...
            '4.0-5-15-2-5-40-5-2-8',...
            '5.0-5-15-2-5-40-5-2-8',...
            '3.0-10-15-2-5-40-5-2-8',...
            '4.0-10-15-2-5-40-5-2-8',...
            '5.0-10-15-2-5-40-5-2-8'};
        classFilePath = '/Volumes/Cumberland/Analysis.noindex/Gauthier/2008-04-30-1/data005-params-union/data005-params-union.classification.txt';
        CELL = 1;
end





% get number of parameter sets
if CELL == 0
    nParams = size(fo,1);
else
    nParams = length(fo);
end

% count number of cells in each param set
count = zeros(nParams,1);
for ii = 1:nParams
    if CELL == 0
        path = [main '/' fo(ii,:) '/' fo(ii,:) '.neurons'];
    else
        path = [main '/' fo{ii} '/' fo{ii} '.neurons'];
    end
    nf = edu.ucsc.neurobiology.vision.io.NeuronFile(path);
    count(ii) = nf.getNumberOfNeurons;
end

% get names of RGC classes (ie on midget, on parasol)
params = struct('cell_type_depth', 1);
cellClasses = load_txt_cell_types(classFilePath, params);

% get sub-classes which group together duplicate copies of each RGC
params = struct('cell_type_depth', 2);
duplicateClasses = load_txt_cell_types(classFilePath, params);

% variable to store all results
results = cell(length(cellClasses),1);

% make variable to identify which cell id comes from which param set
pLookup = [];
for cc=1:nParams
    pLookup = [pLookup repmat(cc,1,count(cc))]; %#ok<AGROW>
end


% in each cell class
for ii=1:length(cellClasses)

    % get class name
    results{ii}.name = cellClasses{ii}.name;


    % identify the index of each unique cell in this class

    % initialize variable storing the indices
    classes = zeros(length(duplicateClasses),1);
    c = 1;

    % go through each duplicate class
    for jj=1:length(duplicateClasses)
        % if the duplicate class contains the name of this cell class
        if ~isempty(strfind(duplicateClasses{jj}.name,cellClasses{ii}.name))
            % note its index
            classes(c) = jj;
            % increment counter
            c = c + 1;
        end
    end
    % c-1 is number of unique cells
    nUnique = c-1;


    % for each unique cell in this class, identify which param sets found
    % it, and how many times

    % create data structures to store results
    results{ii}.cells = zeros(c-1,nParams);
    results{ii}.cell_id = cell(c-1,1);

    % go through each unique cell
    for jj=1:nUnique

        % store the id numbers of each version of this cell
        cid = duplicateClasses{classes(jj)}.cell_ids;
        results{ii}.cell_id{jj} = cid;

        % for each version
        for kk=1:length(cid)

            % identify which param set it belongs to, and
            whichParamSet = pLookup(cid(kk)+1);

            % increment the count for that param set
            results{ii}.cells(jj,whichParamSet) = results{ii}.cells(jj,whichParamSet) + 1;

        end
    end
end
