function [self,cross] = nonoverlapping_metric2(u_spatial_log,mask2)

u_sp_norm = u_spatial_log./repelem(sqrt(sum(u_spatial_log.^2,1)),size(u_spatial_log,1),1);

%% make the shape of filters 2D

sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

metric_su=[];
%figure;
nSU = size(u_sp_norm,2);

u_sp_log = cell(nSU,1);
for isu = 1:nSU
%    subplot(1,size(u_sp_norm,2),isu);
     u_spatial = reshape_vector(u_sp_norm(:,isu),masked_frame,indexedframe);
     u_spatial = u_spatial/max(abs(u_spatial(:)));
     u_sp_log{isu} = u_spatial;
end

%% calculate correlation functions

self=zeros(2*sta_dim1-1,2*sta_dim2-1);
cross = self;
for isu=1:nSU
    for jsu=1:nSU
        xans = xcorr2(u_sp_log{isu},u_sp_log{jsu});

        if(isu==jsu)
            self = self+xans;
            
        else
            cross = cross + xans;
            
        end
        
    end
end
self = self/nSU;
cross = cross/(nSU*(nSU-1));
end