function res=get_ei_neighbors_sorted(electrode,position)
%clockwise, starts at 9 

res=zeros(1,6);

t=find(position(:,1)<position(electrode,1) & position(:,2)==position(electrode,2));
d=(position(t,1)-position(electrode,1)).^2+(position(t,2)-position(electrode,2)).^2;
[~,tt]=min(d);
if ~isempty(tt)
    res(1)=t(tt);
end

t=find(position(:,1)<position(electrode,1) & position(:,2)>position(electrode,2));
d=(position(t,1)-position(electrode,1)).^2+(position(t,2)-position(electrode,2)).^2;
[~,tt]=min(d);
if ~isempty(tt)
    res(2)=t(tt);
end

t=find(position(:,1)>position(electrode,1) & position(:,2)>position(electrode,2));
d=(position(t,1)-position(electrode,1)).^2+(position(t,2)-position(electrode,2)).^2;
[~,tt]=min(d);
if ~isempty(tt)
    res(3)=t(tt);
end

t=find(position(:,1)>position(electrode,1) & position(:,2)==position(electrode,2));
d=(position(t,1)-position(electrode,1)).^2+(position(t,2)-position(electrode,2)).^2;
[~,tt]=min(d);
if ~isempty(tt)
    res(4)=t(tt);
end

t=find(position(:,1)>position(electrode,1) & position(:,2)<position(electrode,2));
d=(position(t,1)-position(electrode,1)).^2+(position(t,2)-position(electrode,2)).^2;
[~,tt]=min(d);
if ~isempty(tt)
    res(5)=t(tt);
end

t=find(position(:,1)<position(electrode,1) & position(:,2)<position(electrode,2));
d=(position(t,1)-position(electrode,1)).^2+(position(t,2)-position(electrode,2)).^2;
[~,tt]=min(d);
if ~isempty(tt)
    res(6)=t(tt);
end

%res=[electrode res];