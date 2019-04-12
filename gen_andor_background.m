%Generate_bacgkround_file
tic
backgroundfile = 'P:\iLocater\AcquisitionCamera\Final\Andor\20190411\Strehl_2_back\Spooled files.sifx'; %backgroundfile

startFrame = 1;

[signal,noFrames,and_size,width,height] = Andor.getDataSetInfo(backgroundfile); %check set exists, frame dimensions, total number of frames)
            
endFrame = noFrames;

dimensions = [height,width];  % assign dimensions to object (used for fitting and plotting)

%%--------------------
% Assign set start/end
%%--------------------

jj = 0;% counter for tracking stored frames

if endFrame > noFrames || endFrame == 0
    endFrame = noFrames;
end

chk1 = round(0.1*(endFrame-startFrame));
chk2 = round(0.5*(endFrame-startFrame));
chk3 = round(0.9*(endFrame-startFrame));


%%-----------------------
% Start processing frames
%%-----------------------

for ii = startFrame:endFrame
    
    if ii == startFrame+chk1
        fprintf('Working on frame %i ...10%% done\n',ii)
    elseif ii == startFrame+chk2
        fprintf('Working on frame %i ...50%% done\n',ii)
    elseif ii == startFrame+chk3
        fprintf('Working on frame %i ...90%% done\n',ii)
    end
    
    %%-----------------------------
    % Read each frame (indiviually)
    %%-----------------------------
    
    %purpose: grabs a single frame and time stamp
    %inputs: flag, size of linearized frame, frame width, frame height, frame number
    %outputs: frame, time stamps
    
    [imageData(:,:,ii),~] = Andor.getFrameInfo(signal,and_size,width,height,ii);
end

    backgroundFrame = median(imageData,3);
    
    save S:\ImageProcessing\CalibrationFiles\Andor\190411_back backgroundFrame
    toc