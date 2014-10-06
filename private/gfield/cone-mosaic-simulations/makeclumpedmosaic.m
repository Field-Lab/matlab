% This file take a set of real LM cone coordinates (x,y only) and makes a
% simulated clumped mosaic, with 3 parameters: pL (average) detlapL (local
% difference in pL,centered around average until either max or min pL
% reaches zero or 1), and scale (FWHM of spatial scale of clumps in pixels)

%type=1 is for L cone, type = -1 is for M cone

%heidi june 2008

function [type actualpL actualdeltapL]=makeclumpedmosaic(LMcoords,pL,deltapL,scale);

numcones=max(size(LMcoords));

%make 512X512 binary white noise image

profile=rand(512,512);

%now filter so that there will be clumps with spatial scale 'scale'

f=fspecial('gaussian',4*scale,scale);
profile=conv2(profile,f,'same');

%now rescale this

av=mean(mean(profile(100:400,100:400)));
maxp=mean(max(profile(100:400,100:400)));
minp=mean(min(profile(100:400,100:400)));

for i=1:512
    for j=1:512
        if profile(i,j)>= maxp
            profile(i,j)=maxp;
        else if profile(i,j) <=minp
                profile(i,j)=minp;
            end
        end
    end
end

range=maxp-minp;

profile=profile-av.*ones(512,512);
profile=av+profile./range; %now we have something where average is about 0.5 and max and min is about 0 and 1

%rescale for close to the desired pL and delta pL (actual values for each
%mosaic will be returned

if deltapL < 1-(2*abs(0.5-pL))
    
    profile=profile-0.5*ones(512,512);
    profile=pL+profile.*deltapL;


else if deltapL >= 1-(2*abs(0.5-pL))
    profile=profile-0.5*ones(512,512);
    profile=pL+profile.*2*(deltapL-0.5+abs(0.5-pL));  
    
    end
   
end


for i=1:512
        for j=1:512
            if profile(i,j)>1
                profile(i,j)=1;
            else if profile (i,j)<0;
                    profile(i,j)=0;
            end
            end
        end
end

% now assign cones based on these local probabilities 1=L -1=M

for i=1:numcones
  if rand(1)> 1-profile(LMcoords(i,1),LMcoords(i,2))
      type(i)=1;
      usedpL(i)=profile(LMcoords(i,1),LMcoords(i,2));
  else
      type(i)=-1;
      usedpL(i)=profile(LMcoords(i,1),LMcoords(i,2));
  end
end

%now get actual average pL and range pL values

actualpL=mean(usedpL);
actualdeltapL=max(max(usedpL))-min(min(usedpL));