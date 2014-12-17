function RGB=RGB_weights(datarun,cell)

% Error if the right part of datarun doesn't exist
try 
    isempty(datarun.vision.timecourses(cell).r);
catch
    error('Requires datarun.vision.timecourses')
end

% make a matrix of the time courses
RGB_matrix=[datarun.vision.timecourses(cell).r ...
    datarun.vision.timecourses(cell).g ...
    datarun.vision.timecourses(cell).b];

% SVD
[U,S,V]=svd(RGB_matrix);

% Take the first component
RGB=V(:,1);

% Warning if the first component doesn't dominate
if S(1,1)/S(2,2) < 10
    warning(['Ratio of first principle component to second is only ' num2str(S(1,1)/S(2,2))])
end

end
