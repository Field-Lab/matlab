j=39:60;
j(15)=[];
basePath='S:\user\alexandra\MEA_data\tmp\MEA_bin\channels\channel';

tmp1=['20121129_CH',int2str(k),'_0307_hlc2_nd%nd_2%_BG1_FWND5ND1.mat'];
tmp=['20121129_CH',int2str(k),'_0309_hlc2_HCseq2_FW0ND2.phys'];
for k=j
    cd([basePath,int2str(k)]);
    a=dir([basePath,int2str(k),'\*.mat']);
    tic
    for i=1:length(a)
        if i==1
            newName=['20121129_CH',int2str(k),'_0309_hlc2_HCseq2_FW0ND2.mat'];
            oldName=a(end).name;
        elseif i==2
             newName=['20121129_CH',int2str(k),'_0307_hlc2_nd%nd_2%_BG1_FWND5ND1.mat'];
             oldName=a(end-1).name;
        else
            newName=a(end-i+3).name;
            oldName=a(end-i+1).name;
        end            
        movefile(oldName,newName);
    end
    toc
end


b=dir([basePath,'37\*mat'])
k=dir([basePath,'38\*mat'])
cd([basePath,'38']);
for i=306:-1:1 
    i
    newName=b(i).name;
    newName=[newName(1:11),'38',newName(14:end)];
    oldName=k(i).name;
    if ~strcmp(oldName,newName)
        movefile(oldName,newName);
    end
end