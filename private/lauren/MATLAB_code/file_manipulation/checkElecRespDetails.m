clear all

%cd /snle/lab/Experiments/Array/Analysis/2010-11-22-5/data008/
cd /Volumes/Ripple/Analysis/lauren/2011-01-11-0/data032/

fileIndex = 0;
files = dir;

%removes non-elecResp files/folders from list
for i = 1:length(files)
    i
    if ~isempty(strfind(files(i).name, 'elecResp')) && isempty(strfind(files(i).name, '._'))
        fileIndex = fileIndex + 1;
        fileNames{fileIndex} = files(i).name; %#ok<SAGROW>
    end
end

clear fileIndex files i

%%

back = text_progress_bar;
%stimVecSizes = zeros(length(fileNames),1);

for xx = 1:length(fileNames)    
    load(fileNames{xx})

    back = text_progress_bar(back,xx/length(fileNames));
    
    %%check details of your choosing
    
    
    stimVecSize1 = size(elecResp.stimInfo.pulseVectors{1});
    for ii = 2:length(elecResp.stimInfo.movieNos)
        if ~all(stimVecSize1 == size(elecResp.stimInfo.pulseVectors{ii}));
            disp(['not all stimulus vectors in ' fileNames{xx} ' are the same size'])
            break
        end
    end
    
    %patternNo = elecResp.stimInfo.patternNo;
    %stimVecSizes(patternNo) = stimVecSize1(1);
    
    

    clear elecResp
end
