function res=robust_mean(data,percent)
% robust_mean 
% res = robust_mean(data,percent)
%

if nargin==1,
    percent=.7;
end

if size(data,1)==1
    data=data';
end

d=sort(data);
nr=ceil(size(d,1)*percent);

if nr<size(d,1)
    for j=1:size(d,2)    
        tem=zeros(1,size(d,1)-nr);
        for i=1:size(d,1)-nr
           tem(i)=d(i+nr,j)-d(i,j); 
        end

        [junk,tmin]=min(tem);

        res(j)=mean(d(tmin:tmin+nr,j));
    end
else
	res=mean(d);
end







if 0
d=sort(data);
nr=ceil(length(d)*percent);

if nr<length(d)
    tem=zeros(1,length(d)-nr);
    for i=1:length(d)-nr
       tem(i)=d(i+nr)-d(i); 
    end

    [junk,tmin]=min(tem);

    res=mean(d(tmin:tmin+nr));
else
    res=mean(d);
end
end