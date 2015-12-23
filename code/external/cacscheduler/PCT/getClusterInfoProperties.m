function props = getClusterInfoProperties()

props = java.util.HashMap();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We pull some properties from the static ClusterInfo class to
    % configure the job.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %The HPC2008 template roughly corresponds to a queue.
    templateName = ClusterInfo.getQueueName();
    if ~isempty(templateName)
        props.put('TEMPLATE_NAME',templateName);
    end
    %This is the estimated walltime for the job
    wallTime = ClusterInfo.getWallTime(true);
    if ~isempty(wallTime)
        props.put('WALLTIME',wallTime);
    end
    %Pull the PPN (This is not currently used!
    ppn = ClusterInfo.getProcsPerNode();
    if ~isempty(ppn)
        props.put('PPN',ppn);
    end
    %Pull the Project Name
    projectName = ClusterInfo.getProjectName();
    if ~isempty(projectName)
        props.put('PROJECT_NAME',projectName);
    end
    %Pull the requested architecture (gpu, bigmem, something)
    arch = ClusterInfo.getArch();
    if ~isempty(arch)
        props.put('ARCHITECTURE',arch);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%