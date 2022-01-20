% Combine LAMMPS trajectory files in MATLAB format to a single .mat file
% E. Weiand - 01/2022

clear
clc

%% parameters

fMin = 1;   % minimum file index
fMax = 2;   % maximum file index


%% processing

sortedDataFull = [];
volDataFull = [];
writeToRestartRatio = 2;
for ff=fMin:fMax
    fileName = ['tabTraj_nvt_pt' num2str(ff) '.mat'];
    
    load(fileName)
    
    len = length(sortedData);
    subtr = mod(len,writeToRestartRatio);
    
    sortedDataFull = [sortedDataFull; sortedData(1:len-subtr)];
    volDataFull = [volDataFull; volData(1:len-subtr)];
end

sortedData = sortedDataFull;
volData = volDataFull;

outName = ['tabTraj_nvt_pt' num2str(fMin) 'to' num2str(fMax) '.mat'];
save(outName,'sortedData','volData','-v7.3')