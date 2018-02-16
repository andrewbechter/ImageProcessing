function [allLocations,data] = Dir_Calc(folderNames,dataStorageLocation)

% N is the number of directories we want to process i.e Folders, Subfolders
% M is the number of files in each directory

pathnames = strcat(dataStorageLocation,folderNames);

N = size(folderNames,2);

for ii = 1:N
    file_id = strcat(pathnames{ii},'\*.tiff');
    allFiles = dir(file_id);
    allNames = {allFiles.name};
    M = size(allNames,2); %M is the number of files in the directory
    for jj= 1:M
        allLocations{jj,ii} = strcat(pathnames{ii},'\',allNames(jj));
        data{jj,ii} = double(imread(char(allLocations{jj,ii})));
    end
end

disp('finished!')

end


