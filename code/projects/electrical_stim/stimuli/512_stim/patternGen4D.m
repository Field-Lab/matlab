function [Array, numPatterns]= patternGen4D(angleIncrement)
%% 4D Pattern Generator

% angleIncrement: given in radians, increment to sample the electrode
%                 stimulus space

% inc = 2*pi/12; 
inc = angleIncrement; 
i = 1; 
theta= 0:inc:(2*pi-inc);
phi= 0:inc:(2*pi-inc);
psi = 0:inc:(2*pi-inc);
% theta= 0:inc:(2*pi-inc);
% phi= 0:inc:(2*pi-inc);
% psi = 0:inc:(2*pi-inc);

numPatterns = length(theta)^3; 
Array = zeros(4,numPatterns); 
for t = 1:length(theta)
    for p = 1:length(phi)
        for ps = 1:length(psi)
            e1 = sin(psi(ps))*sin(phi(p))*cos(theta(t));
            e2 = sin(psi(ps))*sin(phi(p))*sin(theta(t));
            e3 = sin(psi(ps))*cos(phi(p));
            e4 = cos(psi(ps));
            Array(:,i) = [e1 e2 e3 e4];
            i = i+1;
        end
    end
end
% keyboard; 
repeatPer = 2*pi/angleIncrement; 
if repeatPer/2 == round(repeatPer/2)
%     Array(:,repeatPer+1:repeatPer/2:end)
    Array(:,repeatPer+1:repeatPer/2:end)=[];
else
%     Array(:,repeatPer+1:repeatPer:end)
    Array(:,repeatPer+1:repeatPer:end)=[];
end
% k=1;
% while 1 % k = 1:size(Array,2)
% check = Array(:,k); 
% i = 1;
% while 1
%      %:size(Array,2)
%     if all(Array(:,i) == check) && i ~= k
%         Array(:,i) = []; 
%     end
%     if i == size(Array,2)
%         break;
%     end
%     i = i+1;
% end
% if k == size(Array,2)
%     break; 
% end
% k = k+1 
% 
% end
end