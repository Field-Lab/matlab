GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\GaussesFiles\'
Pattern=13;
fid=fopen([GaussesFilesPath 'AllGausses_p' num2str(Pattern)],'r');
        a=fread(fid,'double');
        fclose(fid);    
        b=reshape(a,length(a)/5,5);  