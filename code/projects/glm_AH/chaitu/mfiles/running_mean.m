function m = running_mean(x,sample_length,padval)

    if (~exist('padval','var'))
        padval = 0;
    end

    m = fastconv(x,1/sample_length.*ones(size(x,1),sample_length),size(x,1),size(x,2),padval);
    1;
    %for j=1:sample_length
    %    m(:,j) = mean(x(:,1:j),2);
    %end