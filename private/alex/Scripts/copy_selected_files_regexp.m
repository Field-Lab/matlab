% copy
basedDest='S:\user\alexandra\MEA_data\20120827_1\MEA_bin\channels\';
basedOrig='H:\20120827_1\MEA_bin\channels\';
j=1:60;
j(15)=[];
for k=j
    k
    mainPath=[basedOrig,'channel',int2str(k),'\'];
    origFile=dir(mainPath);
    destPath=[basedDest,'channel',int2str(k),'\'];
    if ~exist(destPath,'dir')
        mkdir(destPath);
    end
    file_list=[];
    for i=1:length(origFile)
        if ~isempty(regexp(origFile(i).name,'ffflicker', 'once'))
            file_list=[file_list i];
        end
    end
    for i=file_list
        copyfile([mainPath,origFile(i).name],[destPath,origFile(i).name]);
    end
end