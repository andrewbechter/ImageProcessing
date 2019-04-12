%Generate_bacgkround_file
tic
backgroundfile = 'P:\iLocater/AcquisitionCamera\Final\FiberChannel\Comb_set6_back\'; %backgroundfile

%backgroundfile = '/Volumes/Projects/iLocater/AcquisitionCamera/Final/Fiber/L2Set4_back_3/';

[~,backNames] = Basler.getDataSetInfo(backgroundfile);

for kk = 1:length(backNames)
    
    [backData] = Basler.getFrame(backgroundfile,backNames(kk));
    
    backgroundFrame(:,:,kk) = backData;
    
end

%sig_back = std(backgroundFrame,1,3);

backgroundFrame = median(backgroundFrame,3);

save CalibrationFiles/14umBasler/comb_set6_back backgroundFrame
toc