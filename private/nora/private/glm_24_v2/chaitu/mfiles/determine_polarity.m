function isup = determine_polarity(spfilt,l)%,m,n)

[idx idx_sorted] = sort(abs(spfilt),'descend');

%[peak_i peak_j] = ind2sub([m n],idx_sorted(1));

%irange = repmat((i-1:i+1)',3);
%jrange = repcols(j-1:j+1,3)';

%idx = sub2ind([m n],irange, jrange);

isup = (mean(spfilt(idx_sorted(1:l))) > 0);