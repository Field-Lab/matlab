classdef (Sealed) ClusterInfo < handle
    %CLUSTERINFO ClusterInfo class.
    %
    %   Provides cluster information about:
    %      QueueName
    %      WallTime
    %
    %   Example:
    %     % Prior to calling a submit script, set a property, as such:
    %  (all equivalent to 10 minute walltime):
    %  ClusterInfo.setWallTime(10); 
    %  ClusterInfo.setWallTime('00:10:00');
    %  ClusterInfo.setWallTime('10m');
    %  ClusterInfo.setWallTime('600s');
    %
    %     % Within the submit script, access the property, as such:
    %     wt = ClusterInfo.getWallTime();   
    
    
    %   Copyright 2009 The MathWorks, Inc.
    %   Raymond S. Norris (raymond.norris@mathworks.com)
    %   $Revision: 71 $ $Date: 2009-06-08 09:28:12 -0400 (Mon, 08 Jun 2009) $
    %   
    %   Modified for TUC
    %   Nate Woody (naw47@cac.cornell.edu)
    %   Cornell Center for Advanced Computing
    %   7/30/2010
    %   Eric Chen (emc256@cornell.edu)
    %   Cornell Center for Advanced Computing
    %   8/13/2010
    
    methods (Static,Access='private')
        function gp = Group()
            gp = 'ClusterInfo';
        end
    end
    
    methods (Access='private')
        function obj = ClusterInfo()
        end
    end
    
    methods (Static)
        function clearClusterInfo()
            try
                rmpref(ClusterInfo.Group)
            catch E %#ok<NASGU>
            end
        end
            
        function setProcsPerNode(ppn)
            if nargin==0 || isnumeric(ppn)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Procs per node must be an integer.')
            end
            setpref(ClusterInfo.Group,'ProcsPerNode',ppn)
        end
        function ppn = getProcsPerNode()
            try
                val = getpref(ClusterInfo.Group,'ProcsPerNode');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = [];
            end
            ppn = val;
        end
        
        
                function setArch(a)
            if nargin==0 || ischar(a)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Arch must be a character string.')
            end
            va = ClusterInfo.ValidArch();
            idx = strmatch(a,va,'exact');
            if isempty(idx)
                error(['Valid arches: ' strtrim(strrep([va{:}],'bit','bit '))])
            end
            ac = ClusterInfo.ArchComplex;
            setpref(ClusterInfo.Group,'Arch',ac{idx})
        end
        function a = getArch()
            try
                val = getpref(ClusterInfo.Group,'Arch');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = [];
            end
            a = val;
        end
        
        
        function setQueueName(qn)
            if nargin==0 || ischar(qn)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Queue name must be a character string.')
            end
            setpref(ClusterInfo.Group,'QueueName',qn)
        end
        function ch = getQueueName()
            try
                val = getpref(ClusterInfo.Group,'QueueName');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = [];
            end
            ch = val;
        end
        
        
        function setProjectName(pn)
            if nargin==0 || ischar(pn)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Project name must be a character string.')
            end
            setpref(ClusterInfo.Group,'ProjectName',pn)
        end
        function pn = getProjectName()
            try
                val = getpref(ClusterInfo.Group,'ProjectName');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = [];
            end
            pn = val;
        end
        function resetWallTime()
            % ClusterInfo.resetWallTime(minutes) - use system default
            % walltime
            rmpref(ClusterInfo.Group,'WallTime')
        end
        function setWallTime(wt)
            % ClusterInfo.setWallTime(minutes) -
            % Prior to calling a submit script, set walltime as such:
            %  (all equivalent to 10 minute walltime):
            %  ClusterInfo.setWallTime(10);
            %  ClusterInfo.setWallTime('00:10:00');
            %  ClusterInfo.setWallTime('10m');
            %  ClusterInfo.setWallTime('600s');
            
            dd = 0;
            hh = 0;
            mm = 0;
            ss = 0;
            if isnumeric(wt)==true
                mm = wt;
            elseif nargin==0
                error('distcomp:clusterinfo:InvalidType', ...
                    'Please specify walltime');
            else
                regex_hourminsec = '(\d\d?\d?):(\d\d?):(\d\d?)';
                regex_name = '(?<day>\d+d)?(?<hour>\d+h)?(?<minute>\d+m)?(?<second>\d+s)?';
                tok = regexp(wt,regex_hourminsec,'tokens');
                
                if ~isempty(tok)
                    dd = 0;
                    hh = str2num(tok{1,1}{1,1});
                    mm = str2num(tok{1,1}{1,2});
                    ss = str2num(tok{1,1}{1,3});
                    % ignore seconds
                    %return
                end
                tok = regexp(wt,regex_name,'names');
                if ~isempty(tok)
                    % convert days/hour/minute/second
                    dd = convertTime(wt,'d');
                    hh = convertTime(wt,'h');
                    mm = convertTime(wt,'m');
                    ss = convertTime(wt,'s');
                    
                    % ignore seconds
                end
                 
            end
            if ss > 60
                mm = mm + ceil(ss/60);
            end
            if mm > 60
                hh = hh + idivide(int32(mm),60,'floor');
                mm = mod(mm,60);
            end
            if hh > 24
                dd = dd + idivide(int32(hh),24,'floor');
                hh = mod(hh,24);
            end
            total_minutes = 24*60*dd + 60*hh + mm;
            if total_minutes==0
                error('invalid walltime specified.  please verify the walltime format.');
            end
            if total_minutes<10
                error('invalid walltime specified.  must be at least 10 minute.');
            end
            fprintf('walltime set to %d day(s), %d hour(s),%d minute(s)\n',dd,hh,mm);
            fprintf('total minutes: %d\n',total_minutes);
            setpref(ClusterInfo.Group,'WallTime',int2str(total_minutes))
            function out = convertTime(str,chr)
                tok = regexpi(str,['(\d+)' , chr],'tokens');
                if isempty(tok)
                    out = 0;
                    return;
                end
                out= str2num(tok{1,1}{1,1});
            end
        end
            
 
        function wt = getWallTime(quiet)
            % ClusterInfo.getWallTime() - returns number of minutes
            % also prints walltime in days/hours/minutes
            try
                val = getpref(ClusterInfo.Group,'WallTime');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = [];
            end
            if isempty(val)
                wt = val;
                if nargin == 0
                    fprintf('walltime set to system default');
                end
                return;
            end
            mm = str2num(val);
            hh = 0;
            dd = 0;
            if mm > 60
                hh = hh + idivide(int32(mm),60,'floor');
                mm = mod(mm,60);
            end
            if hh > 24
                dd = dd + idivide(int32(hh),24,'floor');
                hh = mod(hh,24);
            end
            total_minutes = 24*60*dd + 60*hh + mm;
            if total_minutes==0
                error('invalid walltime specified.  please verify the walltime format.');
            end
            if total_minutes<1
                error('invalid walltime specified.  must be at least 1 minute.');
            end
            if nargin == 0
                fprintf('walltime set to %d day(s), %d hour(s),%d minute(s)\n',dd,hh,mm);
                fprintf('total minutes: %d\n',total_minutes);
            end
            wt = val;
        end
        
        
    end
    
end
