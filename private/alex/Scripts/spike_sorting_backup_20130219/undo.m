function undo(flag,xName,yName)

global rgc ifUpdateISI ifUndo
persistent rgc_backup0 rgc_backup1 rgc_backup2 rgc_backup3 step

switch flag
    case 0
        rgc_backup3=rgc_backup2;
        rgc_backup2=rgc_backup1;
        rgc_backup1=rgc_backup0;
        rgc_backup0=rgc;
        if step>1
            step=step-1;
        end
    case 1
        if step<4
            rgc=rgc_backup1;
            rgc_backup0=rgc;
            rgc_backup1=rgc_backup2;
            rgc_backup2=rgc_backup3;
            rgc_backup3=cell(1,1);
            ifUpdateISI=1;
            ifUndo=1;
            redraw(xName,yName)
            step=step+1;
        else
            display('buffer is empty')
        end
    case 2
        rgc_backup0=cell(1,1);
        rgc_backup1=cell(1,1);
        rgc_backup2=cell(1,1);
        rgc_backup3=cell(1,1);
        step=1;
end

end
