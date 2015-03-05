%% get convolved spike rates from full FFFlicker stimuli

clear
cd('s:\user\alexandra\scripts')
path2save=['S:\data\Anahit\Jimpy\ldMEA\results\Jimpy_New']; % output path
codeWord='FFFlicker';
load('S:\data\Anahit\Jimpy\ldMEA\results\Jimpy_New\forAlex.mat');
NDlist={'8','6','4'}; % NDs to analyze

load(fullfile(path2save,'CellAttributes'))
test=CellAttrib{1,1}{1,1}(:,:);
for o=2:1+numel(NDlist);
   for u=1:3;
       if o==2;
           test3=CellAttrib{1,o}{1,1}(:,:);
           test4=num2cell(CellAttrib{1,o}{1,3}(:,:));
           test5=num2cell(CellAttrib{1,o}{1,2}(:,:));
       elseif o==3;
           test6=CellAttrib{1,o}{1,1}(:,:);
           test7=num2cell(CellAttrib{1,o}{1,3}(:,:));
           test8=num2cell(CellAttrib{1,o}{1,2}(:,:));
       elseif o==4;
           test9=CellAttrib{1,o}{1,1}(:,:);
           test10=num2cell(CellAttrib{1,o}{1,3}(:,:));
           test11=num2cell(CellAttrib{1,o}{1,2}(:,:));
       end
   end
end
forAlex=[test test3 test4 test5 test6 test7 test8 test9 test10 test11];

save(fullfile(path2save,'forAlex'),'forAlex')

for cellID=1:329
    cellID
    pol(cellID,1)=forAlex{cellID,4};
    pol(cellID,2)=forAlex{cellID,7};
    pol(cellID,3)=forAlex{cellID,10};

%     if sum([pol(cellID,:)>=0 pol(cellID,:)<=1])==6 % takes only from cells that have good polarities at all NDs
        date=forAlex{cellID,1};
        if exist(['/mnt/muench_data/data/alexandra/MEA_data/analysis/!WT_Hartwig/',date,'/'],'dir')
            path2load=['/mnt/muench_data/data/alexandra/MEA_data/analysis/!WT_Hartwig/',date,'/'];
             expType(cellID)=1; % wild type
        else            
            path2load=['/mnt/muench_data/data/alexandra/MEA_data/analysis/!KO_Hartwig/',date,'/'];
             expType(cellID)=-1; % KO
        end
        load([path2load,'LF_',codeWord]);
        for k=1:size(names,2)
            if ~isempty(regexp(names{k},forAlex{cellID,2}))
                onOff(cellID)=mean(pol(cellID,:));  
                AllonOff{cellID}=pol(cellID,:);
                load([path2load,'protocols_',codeWord])
                mainpath=['/mnt/muench_data/data/Hartwig/MEA_data/',date,'/'];                
                load([%% selected cells
mainpath,'units/',units(k).name]);
                for i=1:14
                    spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
                    if ~isempty(spikes)
                        tmp=convolved(spikes,40,210000);
                        FiringRate(:,i,cellID)=tmp(121:end-120);
                        tmp=FiringRate(:,i,cellID);
                        for mc=1:10
                            take=correctedProtocols(1200*(mc-1)+1,1,i):correctedProtocols(1200*mc,1,i)-100;
                            tmp1(:,mc)=tmp(take(1:19900));
                        end
                        tmp=reshape(tmp1(:,1:2:end),19900*5,1);
                        full_frSTD_HC(i,cellID)=std(tmp);
                        full_frMean_HC(i,cellID)=mean(tmp);
                        full_fr_HC{i,cellID}=tmp;
                        tmp=reshape(tmp1(:,2:2:end),19900*5,1);
                        full_frSTD_LC(i,cellID)=std(tmp);
                        full_frMean_LC(i,cellID)=mean(tmp);
                        full_fr_LC{i,cellID}=tmp;
                        full_frspont(i,cellID)=mean(FiringRate(250:2000,i,cellID));
                        full_frSTDspont(i,cellID)=std(FiringRate(250:2000,i,cellID));
                        full_fr_spont{i,cellID}=FiringRate(250:2000,i,cellID);
                    end
                end

                break
            end
        end

%     end
    
end
save([path2save,'FR_Full_Range_',codeWord],'full_frMean_HC','full_frSTD_HC','full_frMean_LC','full_frSTD_LC','full_frspont','full_frSTDspont','expType')
save([path2save,'FR_raw_Full_Range_',codeWord],'full_fr_HC','full_fr_LC','full_fr_spont','-v7.3')

