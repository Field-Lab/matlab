function n = howManyClusters(elecs, mainelec)
%%Given active elecs, returns number of clusters of adjacent electrodes

elecs = elecs';
%Load adjacency matrix
load('adj_mat_512.mat');

%seed with mainelec
adjset = adj_mat_512(mainelec); adjset = [adjset{1} mainelec];
elecs = elecs(elecs~=mainelec);
nelecs = length(elecs);

%See if any elecs are in this set
n = 1; %number of clusters
accelecs = 0;
while accelecs < nelecs %while the accounted number of electrodes is less than their total number
	foo = intersect(adjset, elecs); 
	accelecs = accelecs + numel(foo);
	if ~isempty(foo)
		for f = foo
			elecs = elecs(elecs~=f); %remove from list
			grot = adj_mat_512(f);
			adjset = [adjset f grot{1}];
			accelecs = accelecs + 1;
		end
	else
		n = n + 1;
		f = elecs(1);
		elecs = elecs(elecs~=f); %remove from list
		grot = adj_mat_512(f);
		adjset = [adjset f grot{1}];
		accelecs = accelecs + 1;
	end
end
