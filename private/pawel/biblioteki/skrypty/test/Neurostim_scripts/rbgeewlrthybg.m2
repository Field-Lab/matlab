Movie=41;


f=fopen(['D:\Home\Pawel\analysis\2010-09-14-0\data008\mv_' num2str(Movie)]);
        a=fread(f,'int16');
        fclose(f);
        b=reshape(a,NumberOfElectrodes,NumberOfRepetitions,Length);
        s1=reshape(b(1,:,:),NumberOfRepetitions,Length);