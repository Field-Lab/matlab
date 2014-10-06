function add2params(old_file_name,new_file_name,col_names,col_types,data)

%old_file_name='/braid/snle/analysis-archive/Experiments/Array/Analysis/2005-05-02-0/2005-05-02-0/data002/data002-play.params';
%new_file_name='/braid/snle/analysis-archive/Experiments/Array/Analysis/2005-05-02-0/2005-05-02-0/data002/data002-play-new.params';
%old_file_name='/Volumes/Twist/data002/data002-play.params';
%new_file_name='/Volumes/Twist/data002/data002-play-new.params';
%col_names={'test1', 'test2'};
%col_types={'Double', 'DoubleArray'};
%clear data
%data{1,1}=77.77;data{1,2}=88;data{1,3}=99.99;
%data{2,1}=[1.1 2 3 4];data{2,2}=[5 6 7 8.8];data{2,3}=[9 10 11 12.12];
%
%greschner

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile(old_file_name);

columnnames=paramsFile.getColumnNamesAndTypes;
l_columnnames=length(columnnames);
idlist=paramsFile.getIDList;

if ~( size(data,1)==length(col_names) & size(data,1)==length(col_types) & size(data,2)==length(idlist) )
    error('add2params: wrong input size')
end

for i=1:length(col_names)
    columnnames(l_columnnames+i,1)=java.lang.String(col_names{i});
    columnnames(l_columnnames+i,2)=java.lang.String(col_types{i});
end

newparamsFile=edu.ucsc.neurobiology.vision.io.ParametersFile(new_file_name, columnnames, length(idlist));

for i=1:length(idlist)
    r=paramsFile.getRow(idlist(i));

    for ii=1:length(r)
        rr{ii}=r(ii);
    end
    
    for ii=1:length(col_names)
        rr{l_columnnames+ii}=data{ii,i};
    end
    
    newparamsFile.addRow(rr);
    
end


flushsave = true;
newparamsFile.close(flushsave);




