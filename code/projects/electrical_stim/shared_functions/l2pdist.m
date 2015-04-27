function dist = l2pdist(line, point)

x = [line(2,:)-line(1,:);point-line(1,:)];

dist = abs(det(x))/norm(line(2,:)-line(1,:));

end