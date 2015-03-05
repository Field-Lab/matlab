function tracing(src,traceUI)
%
global unitPoly getXdata getYdata axx

if get(traceUI,'Value')==1
    set(traceUI,'String','tracing...','BackgroundColor',[0.2 0.8 0.2]);
    set(src,'WindowButtonDownFcn',@start_trace)
else
    set(src,'WindowButtonDownFcn','')
    set(src,'WindowKeyPressFcn','')
    set(traceUI,'String','press to trace','BackgroundColor',[0.5 0.5 0.5]);
end
unitPoly=[];
acc=[];
ah=gca;
tmp_acc=[];
    function start_trace(src,ev) 
        if strcmp(get(src,'SelectionType'),'normal')
            set(src,'pointer','arrow')
            set(src,'WindowKeyPressFcn',@trace)
        else
            display('SelectionType is wrong')
            return
        end
        
        function trace(src,event)
%             if event.Character == 'x'
%                 axx=[get(gca,'XLim') get(gca,'YLim')];
%                 patch(acc(:,1),acc(:,2),'y','LineWidth',2);
%                 axis(axx)
%                 drawnow
%                 pause(0.2)
%                 unitPoly=[acc(:,1),acc(:,2)];
%                 redraw(get(getXdata,'Value'),get(getYdata,'Value'));
%                 set(src,'WindowButtonMotionFcn','')
%                 set(src,'pointer','arrow')
%             end
            if event.Character == 's'
                tmp_acc=acc;
                axx=[get(gca,'XLim') get(gca,'YLim')];
                for i=1:size(acc,1)-1
                    line(acc([i,i+1],1),acc([i,i+1],2),'color','k','LineWidth',2);
                end
                axis(axx)
                set(src,'pointer','crosshair')
            end
            
            if event.Character == 'c'
                if isempty(tmp_acc)
                    display('Buffer is empty. Press ''s'' first')
                else
                    set(src,'pointer','crosshair')
                    acc=tmp_acc;
                    tmp_acc=[];
                end
            end
        end
        
        cp = get(ah,'CurrentPoint');
        tmp1=[get(gca,'xLim'),get(gca,'yLim')];
        if ~(cp(1,1)>tmp1(1)&&cp(1,1)<tmp1(2)&&cp(1,2)>tmp1(3)&&cp(1,2)<tmp1(4))
            display('Click within plot!')
            cp=[];
            return
        else
            set(src,'pointer','crosshair')
            acc=[cp(1,1) cp(1,2)];
            set(src,'WindowButtonMotionFcn',@mydraw)
        end
        function mydraw(src,evnt)
            cp = get(ah,'CurrentPoint');
            acc=[acc; [cp(1,1) cp(1,2)]];
        end
        set(src,'WindowButtonDownFcn',@di)
        function di(src,~)
            set(src,'pointer','watch')
            axx=[get(gca,'XLim') get(gca,'YLim')];
            patch(acc(:,1),acc(:,2),'y','LineWidth',2);
            axis(axx)
            drawnow
            pause(0.2)
            unitPoly=[acc(:,1),acc(:,2)];
            redraw(get(getXdata,'Value'),get(getYdata,'Value'));
            set(src,'WindowButtonMotionFcn','')
            set(src,'WindowButtonDownFcn',@start_trace)
        end
        
    end

end