classdef gridFTP
%  GRIDFTP is a CAC object that provides simple stateless access to the TUC
%  storage facility to enable you to examine file in your home directory as
%  well as upload and download files to the system.  This tool is primarily
%  useful for uploading or downloading single files or for ensuring the
%  location of files on the storage system.  
%
%  The connection is established and authenticated when the object is
%  instantiated, using X509 certificate.  If you do not currently have a
%  valid credential, you may be prompted to retrieve one from a MyProxy
%  server.
%
%  Note that the ftp connection is stateless, meaning that the connection is
%  always made to your root directory and changing directories is not
%  supported.  Instead, use paths in the list command to examine
%  sub-directories and specify the paths.  
%
%  Directories must be empty in order to be deleted.
%
%  All remotepaths accept either linux or windows path seperators and paths
%  are automatically rooted in your home directory.  You should
%  specify your home using either the Unix '~/' or the windows
%  '\\storage01\matlab\username' paradigm.
%       
%       '~/datadir' == 'datadir' == '\\storage01\matlab\username\datadir'
%       '~'         == ''        == '\\storage01\matlab\username'
% 
% Examples:
%  Create the object and connect to CAC
%       ftp = gridFTP();
%  List the directory contents, always relative to your home
%       ftp.list('');    - list the contents of your home directory
%       ftp.list('MyDataDir');  - list the contents of a subdirectory
%  Download a file to the current directory and delete the remote copy
%       ftp.get('myfile.dat','MyDataDir\resultfile.dat');
%       ftp.delete('MyDataDir\resultfile.dat');
%  Upload a file to a directory that you created in your home directory
%       ftp.createDir('ReferenceData');
%       ftp.put('datafile.txt',''ReferenceData\datafile.txt');
%  Examine the directory to ensure the file was uploaded correctly
%       ftp.list('ReferenceData');
%  Move multiple files with a wildcard.  Destination directories are
%  created if neccessary, wildcards are expected in the source field.
%       ftp.mput('lib\geron*.jar','junk',1);
%       ftp.mget('test','mexjobs/job2/g*');
%
%  See also sendFileToCAC, sendDirToCAC, LittleJohn
%
% Copyright 2010 Cornell Center for Advanced Computing
% 

    properties (SetAccess=protected, Hidden=true )
        ftpClient;
    end
    
    methods (Access=protected)
        function linRootedPath = cleanPath(obj,path)
            %first detect and remote any \\storage01 paths
            if strfind(path,'\\storage01\matlab\')
                %warning('You appear to have specified a \\storage01\matlab\username\somedir path.' + ...
                %    ' We are trying to remove it automatically for you, but next time, just specify' + ...
                %    ' the path without the preceeding components.');
                %they specified the windows path, remote it up to the slash
                %AFTER the username, which means going hunting.
                startinglength = 19;
                path = path(20:end);
                allslashes = strfind(path,'\');
                if numel(allslashes) ~= 0
                    path = path(allslashes(1)+1:end);
                else
                    path = '';
                end
            end
            linPath = strrep(path,'\','/');
            %No entry means root.
            if length(linPath) == 0
                linPath = '~';
            end
            if linPath(1) ~= '~'
                linRootedPath = ['~/' linPath];
            else
                linRootedPath = linPath;
            end
        end
    end
    methods (Access=public)
        function obj = gridFTP()  
            %Create the connnection
            ftp = edu.cornell.cac.tuc.cacscheduler.globus.ftp.SimpleFTPConnection();
            %Insert into a restartable wrapper to restart the connection as
            %it times out.
            obj.ftpClient = edu.cornell.cac.tuc.cacscheduler.globus.ftp.RestartableFTPConnection(ftp);
            obj.connect();
        end
        
        function connect(obj)
            obj.ftpClient.connect();
        end
        
        function bool = isFile(obj,path)
            %FTP.ISFILE('PATH') - Check to see if the request path
            linPath = obj.cleanPath(path);
            if obj.ftpClient.exists(linPath)
                bool =  obj.ftpClient.isFile(linPath);
            else
                bool = 0;
            end
        end
        
        function bool = isDir(obj,path)
            %FTP.ISDIR('PATH') - Check to see if the request path is a
            %directory
            linPath = obj.cleanPath(path);
            if obj.ftpClient.exists(linPath)
                bool = obj.ftpClient.isDirectory(linPath);
            else
                bool = 0;
            end
        end
        
        function createDir(obj,path)
            % FTP.CREATEDIR('DIRNAME') - Create the requested directory on
            % the remote system.  Recursive creation is not supported so all
            % subdirectories must exist before createDir is called.
            %
            %  ftp.createDir('datadir/run001');
            %       will create the run001 folder ONLY if datadir already
            %       exists.
            %
            linPath = obj.cleanPath(path);
            obj.ftpClient.createDirectory(linPath);
        end
        
        function result = list(obj,dir,silent)
            % FTP.LIST('DIRNAME') - DIR/LS on TUC.  Returns a listing of the
            % contents of the specified directory.  Either Unix or Windows
            % style path seperators are supported.
            %  These are equivalent:
            %   ftp.list('datadir/run001');
            %   ftp.list('datadir\run001');
            %
            
           if nargin < 3
               silent=0;
           end
           pth = obj.cleanPath(dir);
           if (isempty(strfind(pth,'*')))
                stuff = obj.ftpClient.infolist(pth);
           else
               %Get the full directory listing, remove filename
               [dir,file,ext] = fileparts(pth);
               pth = obj.cleanPath(dir);
               stuff = obj.ftpClient.infolist(pth);
               if (isempty(stuff))
                   ftp.close();
                   ftp.connect();
                   error('FTP connection appeared to fail, reconnecting!');
               end
               %Filter the list here in matlab
               wc = [file ext];
               %we are likely to have "." in here which we'll need to clean
               reg = strrep(wc,'.','\.');
               reg = strrep(reg,'*','(.+)');
               newstuff = java.util.LinkedList();
               r = stuff.size();
               for i = 0:r-1
                   if regexp(char(stuff.get(i).getFileName()),reg) == 1
                       newstuff.add(stuff.get(i));
                   end
               end
               stuff = newstuff;
           end
           r = stuff.size();
           result = {};
           for i=0:r-1
               lsx = stuff.get(i);
               dstr = char(lsx.get('modify'));
               dstr = datestr(datenum(dstr,'yyyymmddHHMMSS'),'local');
               item.mdate = dstr;
               item.owner = char(lsx.get('unix.owner'));
               item.type = char(lsx.get(lsx.TYPE));
               item.size = str2double(char(lsx.get('size')));
               item.name = char(lsx.getFileName());
               if ~silent
                    fprintf('%4s %6s %8d %8s %s\n', item.type, item.owner,...
                        item.size, item.mdate, item.name);
               end
               result{end+1} = item;
           end
        end
        
        function get (obj,lFile,remoteFile)
            % FTP.GET('LOCALFILE', 'REMOTEFILE')
            % Download the remoteFile to localFile.  
            %
            %   ftp.get(fullfile(pwd,'downfile.mat'),'datadir/result.mat');
            %
            lf = java.io.File(lFile);
            remPath = obj.cleanPath(remoteFile);
            obj.ftpClient.get(lf,remPath);
        end
        
        function mget(obj,lDir, remoteDir, force)
            %MGET(LOCALDIR, REMOTESPEC,FORCE) - multiple get from TUC to the local
            %directory.  
            % localDir is a char specifing the local directory where files
            %should be moved to.  Relative and absolute paths are accepted
            %but relative paths are probably safer.  This directory is
            %created if it doesn't exist and any files are overwritten.
            % remoteDir is a char of the remoteFiles to move.  Wildcard
            % (*) are expected here.  Path should be as used in list.
            % force - set to true to move all files without confirmation.
            %
            % Move all .mat files from TUC in the datadir directory to the
            % local directory test.  Confirm every file move.
            %   ftp.mget('test','datadir/*.mat',0)
            %
            % Move all Task files from the job directory to a local
            % directory. Move without confirming.
            %   ftp.mget('c:\projects\matlab\jobs','Job5\*.mat');
            %
            if nargin < 4
                force = 0;
            end
            %Now we need to find the path referred to in lDir
            [path,name] = fileparts(lDir);
            if isempty(path)
                %name must be a directory in the local directory
                lPath = fullfile(pwd,name);
                mkdir(lPath);
            elseif exist(path,'dir')
                %path must be an absolute path
                lPath = lDir;
                mkdir(lPath);
            else
                error('Can''t figure out the desired destination for the copy, try an absolute path');
            end
            
            %get the listing of the remote Dir
            linPath = obj.cleanPath(remoteDir);
            %We should have something like ~/dir/Task* in linPath, so
            %retrieve the parent dir, fileparts seems to work here
            [path,name,ext] = fileparts(linPath);
            files=obj.list(path,1);
            %Now we want to iterate through each file and see if it passes
            %the WildcardFilter
            lf = java.io.File(fullfile([name,ext]));
            wildcard = edu.cornell.cac.tuc.cacscheduler.globus.ftp.WildCardFilter(lf);
            for i=1:length(files)
                if strcmp(files{i}.type,'file')
                    b = wildcard.accept(lf,files{i}.name);
                    if b
                        %This is a file to move
                        lf = java.io.File(fullfile(lPath,files{i}.name));
                        remPath = [path '/' files{i}.name];
                        %fprintf('get(%s,%s)\n',...
                        %       char(lf.getAbsolutePath),...
                        %       remPath);
                        if force
                            obj.ftpClient.get(lf,remPath);
                        else
                            r = input(sprintf('Move %s?  (y/n):', remPath(3:end)),'s');
                            if strcmpi(r,'y')
                                obj.ftpClient.get(lf,remPath);
                            end
                        end
                    end
                end
            end
        end
            
            
        function delete (obj,remoteFile)
            %FTP.DELETE('REMOTEFILE') - Delete the remote file or
            %directory.  Deleting a non-empty directory is an error and
            %will not succeed.
            %
            %  ftp.delete('datadir/result001.mat')
            %  ftp.delete('datadir');
            %
            remPath = obj.cleanPath(remoteFile);
            if obj.ftpClient.isDirectory(remPath)
                obj.ftpClient.deleteDir(remPath);
            else
                obj.ftpClient.deleteFile(remPath);
            end
        end
        
        function deleteeverything(obj,remoteDir,verbose)
            %FTP.DELETEEVERYTHING('remoteDir') - Recursively delete the
            %remote directory (rm -rf).  Wildcards are not accepted.
            % if verbose == true, each deleted file and directory is
            % printed to standardout.
            %
            %Delete all files and subdirectories of a directory, then delete  
            %the directory:
            % ftp.deleteeverything('datadir'); 
            % ftp.deleteverything('somefile.mat');  %Error, use ftp.delete
            if nargin < 3
                verbose = 0;
            end
            remPath = obj.cleanPath(remoteDir);
            if ~obj.ftpClient.exists(remPath)
                error('Remote directory doesn''t appear to exist.');
            end
            if ~obj.ftpClient.isDirectory(remPath)
                error('Remote directory appears to be a file, use FTP.DELETE');
            end
            files = obj.ftpClient.infolist(remPath);
            r = files.size();
            for i=0:r-1
               lsx = files.get(i);
               type = char(lsx.get(lsx.TYPE));
               %If type is a dir, we just call deleteEverything with the
               %new path
               if strcmp(type,'dir')
                   obj.deleteeverything([remPath '/' char(lsx.getFileName())],verbose);
               elseif strcmp(type,'file')
                   if verbose
                       fprintf('deleting %s\n',[remPath '/' char(lsx.getFileName())]);
                   end
                   obj.delete([remPath '/' char(lsx.getFileName())]);
               end
            end
            %Lastly, delete this directory
            if verbose
                fprintf('deleting %s\n',remPath);
            end
            if strcmp(remPath,'~')
                warning('Can''t delete home directory, skipping...');
            else
                obj.delete(remPath);
            end
        end
        function close(obj)
            obj.ftpClient.close();
        end
        function put(obj,lFile,remoteFile)
            %FTP.PUT('LOCALFILE','REMOTEFILE') - Upload the specified file
            %to TUC.  
            %
            %   ftp.put('datafile.txt','remoteDir/thefile.mat');
            %
            remPath = obj.cleanPath(remoteFile);
            lf = java.io.File(lFile);
            obj.ftpClient.put(lf,remPath);
        end
        
        function mput(obj,lFile,remoteDir,force)
            %FTP.MPUT('DATA\*.M','TUCDATADIR', FORCE) - Upload multiple files to
            %TUC.  If force is true then the user is not prompted before
            %every matching file and all matching files are moved.
            % 
            if nargin < 4
                force = false;
            end
            %lf is likely a wildcard type thing with a relative path, sort
            %out an absolute path from it.
            [path,name,ext] = fileparts(lFile);
            %Java.io.File has a different view of relative paths than
            %matlab, so we have to root this based on matlabs pwd
            tPath = fullfile(pwd, path);
            %See if tPath exists,
            if exist(tPath)~= 7
                %Maybe this was an absolute path to begin with?
                tPath2 = fullfile(path);
                if exist(tPath2) ~= 7   
                    error('Can''t resolve the requested local directory.  Tried %s and %s\n', tPath, tPath2);
                end
                tPath = tPath2;
            end
            %Now add the wildcard back in to generate the wildcard string
            lf = java.io.File(fullfile(tPath,[name,ext]));
            wildcard = edu.cornell.cac.tuc.cacscheduler.globus.ftp.WildCardFilter(lf);
            localDir = lf.getParentFile();
            files = localDir.list(wildcard);
            [numFiles,dumb] = files.size();
            fprintf('Found %d matching files in %s\n', numFiles,char(localDir.getCanonicalPath() ));
            %Now create the remote directory if needed
            linPath = obj.cleanPath(remoteDir);
            if obj.ftpClient.exists(linPath)
                if obj.ftpClient.isDirectory(linPath)
                    fprintf('Remote directory %s already exists, using this directory.\n' , linPath);
                else obj.ftpClient.isFile(linPath)
                    error('Requested directory %s is a file, can''t use as a destination.', linPath);
                end
            else
                success = obj.ftpClient.createDirectory(linPath);
                if success
                    fprintf('Remote directory %s created.\n', linPath);
                else
                    error('Unable to create requested remote directory %s', linPath);
                end
            end
            %Now actually move the files
            for i = 1:files.size()
                source = fullfile(char(localDir),char(files(i)));
                dest = obj.cleanPath([linPath '/' char(files(i))]);
                if ~force
                    r = input(sprintf('Move %s?  (y/n):', source),'s');
                    if strcmpi(r,'y')
                        obj.ftpClient.put(java.io.File(source),dest);
                    end
                else
                    obj.ftpClient.put(java.io.File(source),dest);
                end
                %fprintf('ftpClient.put(%s,%s);\n',source,dest);
            end
        end
    end
end