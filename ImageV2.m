classdef ImageV2
    properties
        path
        type
        totalFrames
        Frames
        fileID
    end
    methods
        
        
        function [obj] = ImageV2(path,type)
            
            %tells us what kind of data we are using
            %path to data for processing
            %figures out how many files there are
            % if there are normal amount it reads in the data set.
            % if we need to run in low memory then we have a special
            % function
            
            if nargin < 1
                disp('you must specify a path')
                return
                
            elseif nargin < 2 || strcmp(type,'.tiff')==1
                type = '.tiff'; % need to change this to a check
                obj.fileID = strcat(path,'\*',type);
                allFiles = dir(obj.fileID);
                allNames = {allFiles.name};
                obj.path = path;
                obj.totalFrames = size(allNames,2); %M is the number of files in the directory
                obj = readData(obj);
            elseif strcmp(type,'.sifx')==1
                %do andor files
                %add spooled files to path
            elseif strcmp(type,'.fits')==1 || strcmp(type,'.fit')==1
                %do fits read
                %path is a directory name
                
            end
            
        end
        
        function [obj] = readData(obj)
            allData = cell(obj.totalFrames,1);
            allFiles = dir(obj.fileID);
            allNames = {allFiles.name};
            for jj= 1:obj.totalFrames
                allData{jj,1} = strcat(obj.path,'\',allNames(jj));
                data = double(imread(char(allData{jj,1})));
                obj.Frames(:,:,jj) = data;
            end
        end
 
    end
end
