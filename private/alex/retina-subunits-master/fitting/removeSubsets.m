function I_sc = removeSubsets(I_sc)

% given a set of assignments from cones to subunits,
% check whether any subunit is a strict subset of another
% and remove it. does it recursively, so we end up with a 
% full set of subunits with no repeats

for is_1 = 1:size(I_sc,1)
    for is_2 = is_1+1:size(I_sc,1)
        if sum(I_sc(is_1,:)&I_sc(is_2,:)) == min(sum(I_sc(is_1,:)),sum(I_sc(is_2,:)));
            inds = ones(size(I_sc,1),1);
            if sum(I_sc(is_1,:)) < sum(I_sc(is_2,:))
                inds(is_1) = 0;
                I_sc = I_sc(logical(inds),:);
            else
                inds(is_2) = 0;
                I_sc = I_sc(logical(inds),:);
            end     
            I_sc = removeSubsets(I_sc);
            break
        end
    end
end
                