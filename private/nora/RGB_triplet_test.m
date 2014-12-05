cell=2;

hold on
plot(datarun.vision.timecourses(cell).r,'r')
plot(datarun.vision.timecourses(cell).g,'g')
plot(datarun.vision.timecourses(cell).b,'b')

RGB_linear=[1, linsolve(datarun.vision.timecourses(cell).r,datarun.vision.timecourses(cell).g),...
linsolve(datarun.vision.timecourses(cell).r,datarun.vision.timecourses(cell).b)];

%%
hold on
scatter(datarun.vision.timecourses(cell).r,datarun.vision.timecourses(cell).g,'g')
scatter(datarun.vision.timecourses(cell).r,datarun.vision.timecourses(cell).b,'b')
vec=[min(datarun.vision.timecourses(cell).r) max(datarun.vision.timecourses(cell).r)];
plot(vec, RGB_linear(2)*vec,'g');
plot(vec, RGB_linear(3)*vec,'b');

%% 
RGB=[datarun.vision.timecourses(cell).r ...
    datarun.vision.timecourses(cell).g ...
    datarun.vision.timecourses(cell).b];

[U,S,V]=svd(RGB);
RGB_SVD=V(:,1);