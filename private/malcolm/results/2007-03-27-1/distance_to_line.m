onp=(double(results(:,2))>11 & double(results(:,5))==1);
onm=(double(results(:,2))>11 & double(results(:,5))==3);
dataset1 = [log(double(results(onp,7))),log(double(results(onm,7)))];

onp=(double(results(:,2))<=11 & double(results(:,5))==1);
onm=(double(results(:,2))<=11 & double(results(:,5))==3);
dataset2 = [log(double(results(onp,7))),log(double(results(onm,7)))];

n = repmat([1/sqrt(2),1/sqrt(2)],length(dataset1),1);
distances1 = dataset1-repmat(sum(dataset1.*n,2)*(1/sqrt(2)),1,2);
distances1 = sqrt(sum(distances1.^2,2));
n = repmat([1/sqrt(2),1/sqrt(2)],length(dataset2),1);
distances2 = dataset2-repmat(sum(dataset2.*n,2)*(1/sqrt(2)),1,2);
distances2 = sqrt(sum(distances2.^2,2));