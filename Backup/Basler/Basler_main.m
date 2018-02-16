% Basler Main Script
clear all
close all

% dataStorageLocation = 'P:\iLocater\iLocater_Demonstrator\Lab_Stability\Set22\';
dataStorageLocation = 'P:\iLocater\QuadCell\';

%Takes a user input and displays all the folders found in that directory
%an example is given to show how to only pick certain folder names (i.e
%with the letters MIC in the name)

%%The script will seach in all subfolders of the directory you choose.
%%Files will still be grouped by subfolder in data cell array

if isempty(dataStorageLocation)
    disp('this is not the directory you are looking for...')
    return
else
    %change this line to adjust how the subfolders are chosen.
    %*. will select all subfolders in the directory but you can change that to match a specific identifier
    file_tag = '\*.';      
    datafolders = dir(strcat(dataStorageLocation,file_tag));
    datafolders(1)=[]; %removes hidden '.' folder from the list
    datafolders(1)=[]; %removes hidden '..' folder from the list
    folderNames = {datafolders.name}; % convert to cell array froms structure type
    numsubfolders = size(folderNames,2);% number of rows(subfolders) in directory
    
    [allLocations,data] = Dir_Calc(folderNames,dataStorageLocation);
   
     p = Centroid(data);
     
     h = cell2mat(p(:,:,1));
     v = cell2mat(p(:,:,2));
     
     figure
     hold on
     axis equal
     for ii = 1:size(v,2)
         plot(h(:,ii),v(:,ii),'.')
         stat_v(1,ii) = mean(v(:,ii));
         stat_h(1,ii) = mean(h(:,ii));
         stat_v(2,ii) = std(v(:,ii));
         stat_h(2,ii) = std(h(:,ii));
         stat_v(3,ii) = max(v(:,ii))-mean(v(:,ii));
         stat_h(3,ii) = max(h(:,ii))-mean(h(:,ii));
         
     end
     
%      max_disp_v = max(disp_v);% finds peak displacement from mean
%      max_disp_h = max(disp_h);
     
end