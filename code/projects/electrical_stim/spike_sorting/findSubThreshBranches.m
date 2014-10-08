function branches = findSubThreshBranches(links, thresh, nBranches)

% finds the largest branches (in terms of number of leaves) with link size below the given threshold,
% and returns the IDs of the leaves within them
%
% usage: branches = findSubThreshBranches(links, thresh, nBranches)
%
%
% arguments:
%   links = output of matlab's linkage function
%   thresh = threshold link length
%   nBranches = number of largest branches to find (default = 2)
%
% returns:
%   branches{i} = leaf IDs of ith-largest branch (may be empty if no links are below threshold)


if nargin < 3
    nBranches = 2;
end

% determining largest single branch after each new link
nLink = size(links, 1);


branchSize = zeros(1, 2*nLink + 1);
branchSize(1:nLink+1) = ones(1, nLink+1);

largestBranchID = zeros(1, nLink);

leafIDInBranch = cell(1, 2*nLink+1);
for k = 1:nLink
    leafIDInBranch{k} = k;
end

for k = 1:nLink
    branchSize(nLink+1+k) = branchSize(links(k,1))+branchSize(links(k,2));
    largestBranchID(k) = find(branchSize==max(branchSize),1);
    leafIDInBranch{nLink+1+k} = [leafIDInBranch{links(k,1)} leafIDInBranch{links(k,2)}];
end



mainBranchID = zeros(nBranches, 1);
mainBranchSize = zeros(nBranches, 1);
branches = cell(nBranches, 1);

threshLink = find(squeeze(links(:,3)) > thresh, 1) - 1;

if threshLink > 0
    mainBranchID(1) = largestBranchID(threshLink);
    mainBranchSize(1) = length(leafIDInBranch{mainBranchID(1)});
    branchesToCheck = ones(2*nLink + 1, 1);
    branchesToCheck([links(mainBranchID(1)-(nLink+1), 1:2) mainBranchID(1)]) = 0;
    
    for i = 2:nBranches
    
        for k = threshLink + nLink + 1 : -1 : 2 + nLink %all branches below threshold
            if branchesToCheck(k) %check if current branch is distinct from other identified main branches
                nLeavesInBranch = length(leafIDInBranch{k});
                if nLeavesInBranch > mainBranchSize(i)
                    mainBranchID(i) = k;
                    mainBranchSize(i) = nLeavesInBranch;
                end
            else
                branchesToCheck(links(k-(nLink+1), 1:2)) = 0;
            end
        end

        if mainBranchID(i) ~= 0
            branches{i} = leafIDInBranch{mainBranchID(i)};
            branchesToCheck([links(mainBranchID(i)-(nLink+1), 1:2) mainBranchID(i)]) = 0;
        else %no more distinct branches below threshold
            for j = i:nBranches
                branches{j} = [];
            end
            break
        end
    end
    branches{1} = leafIDInBranch{mainBranchID(1)};
elseif all(squeeze(links(:,3)) < thresh) %all branches below threshold
    branches{1} = 1:(nLink+1);
else % no branches below threshold
    branches{1} = [];
end


















