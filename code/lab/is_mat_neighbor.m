function res=is_mat_neighbor(mat,start, varargin)
%test recursively for contiguity
%
% eg: [t1 t2]=max(mat(:)); res=is_mat_neighbor(mat>t1*thr,t2);
%
%greschner


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('number_of_neighbors', 1);% 
    p.addParamValue('diagonals', 1);% 
    p.parse(varargin{:});
    params = p.Results;


if length(start)==1
    [start(1) start(2)]=ind2sub(size(mat),start);
end

res=zeros(size(mat));
res(start(1),start(2))=1;


res=rev(mat,start,start,res,params.number_of_neighbors,params.diagonals);

if params.number_of_neighbors==2
    res(res<0)=0;
    tres=zeros(size(mat));
    tres(start(1),start(2))=1;
    res=rev(res,start,start,tres,1,params.diagonals);
end

res=logical(res);




% Helper function
function res=rev(mat,actual,start,res,nr,diagonals)

    if diagonals
        t=[1 -1;1 0;1 1;0 -1;0 1;-1 -1;-1 0;-1 1];
    else
        t=[0 -1;0 1;-1 0;1 0];
    end

    for i=1:length(t)
        if (actual(1)+t(i,1))>=1 && (actual(1)+t(i,1))<=size(mat,1) && (actual(2)+t(i,2))>=1 && (actual(2)+t(i,2))<=size(mat,2)
            if mat(actual(1)+t(i,1),actual(2)+t(i,2)) && ~res(actual(1)+t(i,1),actual(2)+t(i,2))
                res(actual(1)+t(i,1),actual(2)+t(i,2))=1;
                res=rev(mat,[actual(1)+t(i,1),actual(2)+t(i,2)],start,res,nr,diagonals);
                if nr==2 && ~all(actual==start)
                    tmat=res;
                    tmat(actual(1),actual(2))=0;
                    tres=zeros(size(tmat));
                    tres(actual(1),actual(2))=1;
                    ttres=conection_(tmat,[actual(1)+t(i,1),actual(2)+t(i,2)],start,tres,diagonals);

                    if ~all(ttres(:))
                        res(actual(1)+t(i,1),actual(2)+t(i,2))=-1;
                    end
                end
            end
        end
    end





           

% Helper function
function res=conection_(mat,start,target,res,diagonals)

    if all(start==target)
        res(:,:)=1;
    end

    if diagonals
        t=[1 -1;1 0;1 1;0 -1;0 1;-1 -1;-1 0;-1 1];
    else
        t=[0 -1;0 1;-1 0;1 0];
    end

    for i=1:length(t)
        if (start(1)+t(i,1))>=1 && (start(1)+t(i,1))<=size(mat,1) && (start(2)+t(i,2))>=1 && (start(2)+t(i,2))<=size(mat,2)
            if mat(start(1)+t(i,1),start(2)+t(i,2)) && ~res(start(1)+t(i,1),start(2)+t(i,2))
                res(start(1)+t(i,1),start(2)+t(i,2))=1;
                res=conection_(mat,[start(1)+t(i,1),start(2)+t(i,2)],target,res,diagonals);
            end
        end
    end


























