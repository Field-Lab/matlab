function mov = activity_movie(datarun, output_filename, varargin)
% activity_movie     spike activity on electrodes e.g wave observation
%
% arguments:  datarun - datarun struct 
%             output_filename
%
%             params    array_type - 61 | 512 | 519,  default 512
%                       start - default 0 sec
%                       stop - duration neurons file
%                       bin - default 1 sec
%                       colormap - default jet
%                       saturate - collapse top % of color values, default 0
%                       log_color - log scale for color, default 0 
%                       marker_size - default 24
%                       fps - frames per sec, playing speed of final movie, default 10fps 
%
%             output    movie
%
% greschner


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('array_type', 512);% 
    p.addParamValue('start', 0);% 
    p.addParamValue('stop', datarun.duration);%
    p.addParamValue('bin', 1);%
    p.addParamValue('colormap', 'jet');%
    p.addParamValue('saturate', 0);%
    p.addParamValue('log_color', 0);%
    p.addParamValue('marker_size', 24);%
    p.addParamValue('fps', 10);%
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;
    
    disp('activity movie - running');
    
    datarun=load_electrode_position(datarun, 'array_type', params.array_type);
    
    %combine spikes from channels
    channels{size(datarun.position,1)}=[];
    for i=1:length(datarun.cell_ids)
        channels{datarun.channels(i)}=[channels{datarun.channels(i)}; datarun.spikes{i}];
    end
    
    %calculate histograms 
    mat=zeros(ceil((params.stop-params.start)/params.bin),size(datarun.position,1));
    for i=1:size(datarun.position,1)
        if ~isempty(channels{i})
            mat(:,i)=histc(channels{i},[params.start:params.bin:params.stop-params.bin]);
        end
    end

    %normalize for color
    col=colormap(params.colormap);
    if params.log_color
        mat=log(mat+1);
    end
    if isempty(max(mat(:)))
        error('no spikes in selected interval');
    end
    mat=mat/max(mat(:))/(1-params.saturate)*63;
    mat=round(mat);
    mat(find(mat>63))=63;
    mat=mat+1;

    %plot frames
    h=figure('MenuBar','none', 'Color', [1 1 1],'Visible','off');
    limits=[min(datarun.position(:,1)) max(datarun.position(:,1));min(datarun.position(:,2)) max(datarun.position(:,2))]*1.1;
    clear mov
    for i=1:size(mat,1)
        clf
        set(gca,'XTick',[],'YTick',[],'XLim',limits(1,:),'YLim',limits(2,:));
        axis equal;
        axis off; 
        hold on
        for ii=1:size(datarun.position,1)
            plot(datarun.position(ii,1),datarun.position(ii,2),'.','color',col(mat(i,ii),:),'MarkerSize',params.marker_size); 
        end
        title(sprintf('%0.2f sec',(i-1)*params.bin+params.start));
        mov(i)=getframe(h);
    end
    
    disp(' - saving movie');
    movie2avi(mov,output_filename,'FPS',params.fps);

end


