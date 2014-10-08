function slider = make_loop_slider_list(start_index,index_min,index_max)
% GUI slider for loop over plot for manual inspection
%
%    start_index=1; index_min=1; index_max=10;
%    mat=rand(100,index_max);
%    
%    sliderFig = figure;
%    slider = make_loop_slider_list(start_index,index_min,index_max);
%    while ishandle(sliderFig)
%       k = round(get(slider,'Value'));
%       plot(mat(:,k));
%       uiwait;
%    end
%
% 2009 greschner
% 2009-06  gauthier, simplified for no buttons
%


slider = uicontrol(gcf,...
    'Style','slider',...
    'Min' ,index_min,'Max',index_max, ...
    'Units','normalized', ...
    'Position',[0,0,0.96,.04], ...
    'Value', start_index,...
    'SliderStep',1/(index_max-index_min) * [1 1],...
    'CallBack', 'uiresume;');

set(gcf,'Toolbar','figure')
