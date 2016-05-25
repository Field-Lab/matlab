function [ dist ] = PD_LineDist( v1, edge1, v2, edge2 )
    
    d = zeros(2,2);
    v = [v1 v2];
    edge = zeros(2,2,2);
    edge(:,:,1) = edge1';
    edge(:,:,2) = edge2';
    
    for e = 1:2 % loop for edges
        for line = 1:2 % loop for line
            p1 = edge(:,floor((line+1)/line),e) - edge(:,line,e);
            p2 = v(:,line)*(v(:,line)'*p1);
            d(e,line) = PD_Dist(p1,p2);
        end
    end
   
    dist = max(max(d));
end

