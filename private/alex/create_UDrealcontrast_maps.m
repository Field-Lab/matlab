savedMap=dlmread('/Volumes/Data/2012-08-21-0/Visual/freeman4_d01_localmax/map-0000.txt');
figure
imagesc(savedMap)

cone_weights=[1 0.8 0.5 0.5];
contrast_list=[];
for i=1:3
    for j=2:4
        if i~=j
            defseq=[0 0 0 0];
            minWeight=min(cone_weights([i j]));
            defseq([i j])=1./cone_weights([i j])*minWeight;
            defseq(i)=-defseq(i);
            contrast_list=[contrast_list; defseq];
            defseq([i j])=-defseq([i j]);
            contrast_list=[contrast_list; defseq];
        end
    end
end

nTrials=30;
contrast_list=repmat(contrast_list,nTrials,1);
tmp=randperm(size(contrast_list,1));
contrast_list=contrast_list(tmp,:);

fid=fopen('/Users/alexth/Desktop/test/tmp.txt','w');
fprintf(fid,'(:n 4 :map "cancellation_stimulus:test_control:map-0000.txt")\r');
for i=1:size(contrast_list,1)-1
    fprintf(fid,'(%f\t%f\t%f\t%f)\r\n',contrast_list(i,:));
end
fprintf(fid,'(%f\t%f\t%f\t%f)',contrast_list(end,:));
fclose(fid)
