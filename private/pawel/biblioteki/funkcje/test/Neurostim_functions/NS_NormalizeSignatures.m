function TracesNew=NS_NormalizeSignatures(Traces,range);

TracesNew=Traces;
ss=size(Traces);
for i=1:ss(1)
    %signature=reshape(Traces(i,:,),ss(2),ss(3));
    for j=1:ss(2)
        data=reshape(Traces(i,j,:),1,ss(3))
        size(data)
        d=abs(data)
        size(d)
        max(max(d))
        %data=data-mean(data);
        TracesNew(i,j,:)=data/max(abs(data))*range;
    end
end
    