function input = tempFilterCones(input,timeCourse,mode)

% first dimension should be duration, second stimulus index
% timeCourse just needs to be a vector
%%
if ~exist('mode','var') || isempty(mode)
    mode = 'full';
end

dur = size(input,1);
nStim = size(input,2);
filt = flipud(timeCourse(:));
for i=1:nStim  
    switch mode
        case 'full'
            filtInput = conv(input(:,i),filt,'full');
        case 'circular'
            filtInput = cconv(input(:,i),filt,length(1:dur));
    end
    
    input(:,i) = filtInput(1:dur);
end