%freqextractiontest
clear 

% function [power,seconds] = CPLPLPowerData(filename)

%%%%%%%----------------TEST file location-----------------%%%%%
% filename = 'S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\pm_04-17-2017_10h54m38s.txt';
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Demonstrator\2016_04_18\72Leo\YbandPower\pm_04-17-2016_22h31m49s.txt';
%%%-----------------File read in---------------%%%%%
format long
fid = fopen(strcat(filename));
tot_num_lines = 0;
title_line_skip = 0;
tline = fgetl(fid);
used_lines = 0;
num_lines = 0;

while ischar(tline)
    if(isempty(strfind(lower(tline),lower('Host name'))) == 0)
        disp('Skipping line')
    elseif(isempty(strfind(lower(tline),lower('Device:'))) == 0)
        disp('Skipping line')
    elseif(isempty(strfind(lower(tline),lower('Date:'))) == 0)
        disp('Skipping line')
    elseif(isempty(strfind(lower(tline),lower('Opened the file at'))) == 0)
        disp('Skipping line')
    elseif(isempty(strfind(lower(tline),lower('Started'))) == 0)
        disp('Skipping line')
    elseif(isempty(strfind(lower(tline),lower('W'))) == 0)
        %disp('Skipping line')
    elseif(isempty(tline) == 1)
        disp('Blank Line')
    elseif(isempty(strfind(lower(tline),lower('Active Channel'))) == 0)
        disp('Skipping line')
    elseif(isempty(strfind(lower(tline),lower('Finished taking the data at'))) == 0)
        disp('Skipping line')
    elseif(isempty(strfind(lower(tline),lower('Program'))) == 0)
        disp('Skipping line')
    else
        tot_num_lines = tot_num_lines + 1;
    end
    
    tline = fgetl(fid);
 
end

frewind(fid)
tline = fgetl(fid);


while ischar(tline)
    
    if(isempty(strfind(lower(tline),lower('Host name'))) == 0)
        disp('Skipping line')
        title_line_skip = title_line_skip + 1;
    elseif(isempty(strfind(lower(tline),lower('Device:'))) == 0)
        disp('Skipping line')
        title_line_skip = title_line_skip + 1;
    elseif(isempty(strfind(lower(tline),lower('Date:'))) == 0)
        disp('Skipping line')
        title_line_skip = title_line_skip + 1;
    elseif(isempty(strfind(lower(tline),lower('Opened the file at'))) == 0)
        disp('Skipping line')
        title_line_skip = title_line_skip + 1;
    elseif(isempty(strfind(lower(tline),lower('Started'))) == 0)
        disp('Skipping line')
        title_line_skip = title_line_skip + 1;
    elseif(isempty(strfind(lower(tline),lower('Active Channel'))) == 0)
        disp('Skipping line')
        title_line_skip = title_line_skip + 1;
    elseif(isempty(strfind(lower(tline),lower('W'))) == 0)
        title_line_skip = title_line_skip + 1;
    elseif(isempty(tline) == 1)
        disp('Blank Line')
        title_line_skip = title_line_skip + 1;
    elseif(isempty(strfind(lower(tline),lower('Finished taking the data at'))) == 0)
        disp('Skipping line')
        title_line_skip = title_line_skip + 1;
    elseif(isempty(strfind(lower(tline),lower('Program'))) == 0)
        disp('Skipping line')
        title_line_skip = title_line_skip + 1;

    else
            
        raw_data = strsplit(tline);
        if(size(raw_data,2) == 3)
            used_lines = used_lines + 1; 
            time_cat=strcat(raw_data(1,1),'-',raw_data(1,2));
            date = datevec(time_cat,'mm/dd/yyyy-HH:MM:ss.FFF');
            seconds(used_lines,1) = datenum(date)*86400;
            power_raw(used_lines,1)=raw_data(1,3);
        else
            disp('Skipping line')
        end        
       
    end
    
    tline = fgetl(fid);
end
exp_time =seconds(used_lines,1)-seconds(1,1);
t(:,1) = seconds-(seconds(1,1));
power=str2num(cell2mat(power_raw));


% dateFormat = 0;
% datetick('x',dateFormat,'keepticks')
%%fitting input or output data to a straight line
%x = seconds;
% y = power;
% p = polyfit(x,y,1)
% yfit = polyval(p,x);
% yfit =  p(1) * x + p(2);
% yresid = y - yfit;
% SSresid = sum(yresid.^2);
% SStotal = (length(y)-1) * var(y);
% rsq = 1 - SSresid/SStotal


%end

   
        
        
        
