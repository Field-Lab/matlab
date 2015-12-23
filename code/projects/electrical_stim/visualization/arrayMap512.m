function arrayMap512(positions)

if nargin < 1
    % Load matrix containing the electrode numbers for the 512-electrode MEA
    temp = load('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/freddy/512elecpositions.mat'); % Find a more general location for this or call a different text file.
    positions = temp.positions;
end
    
figure('position', [0 500 900 500]); scatter(positions(:,1),positions(:,2),200,[0.5,0.5,1],'filled'); axis image; whitebg('black'); axis off; set(gcf, 'InvertHardCopy', 'off');
for x = 1:512
    text(positions(x,1),positions(x,2),num2str(x),'HorizontalAlignment','center');
end

end