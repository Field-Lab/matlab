function arrayMap512(positions)

if nargin < 1
    % Load matrix containing the electrode numbers for the 512-electrode MEA
    temp = load(fullfile(matlab_code_path,...
        'code/projects/electrical_stim/resources/arrayPositions512.mat')); 
    positions = temp.positions;
end
    
figure('position', [0 500 900 500]); 
scatter(positions(:,1),positions(:,2),200,[0.5,0.5,1],'filled'); 
axis image; axis off; 
for x = 1:512
    text(positions(x,1),positions(x,2),num2str(x),'HorizontalAlignment','center');
end

end