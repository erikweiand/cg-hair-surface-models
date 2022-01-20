% Read trajectory files from LAMMPS for post-processing
% E. Weiand - 01/2022

clear
clc

% minimum and maximum trajectory part files
nMin=1;
nMax=1;

for ff=nMin:nMax
    % read file
    fileName = ['traj_nvt_pt' num2str(ff) '.lammpstrj'];
    filetext = fileread(fileName);

    % split by timestep delimiter
    splitText = strsplit(filetext, 'ITEM: TIMESTEP');
    splitText = splitText(2:end)';

    len = length(splitText);
    sortedData = cell(len,1);
    volData = cell(len,1);

    for i=1:len
        % split by lines
        currStr = splitText{i};
        currStr = splitlines(currStr);
        currVolStr = currStr(6:8);
        currStr = currStr(10:end-1);

        % volume data
        currVolData = cellfun(@(x) strsplit(x, ' '), currVolStr, 'UniformOutput', false);
        currVolData = vertcat(currVolData{:});
        currVolData = cellfun(@str2num, currVolData);
        volData{i} = currVolData;

        % split lines by delimiter and convert to num matrix
        currData = cellfun(@(x) strsplit(x, ' '), currStr, 'UniformOutput', false);
        currData = vertcat(currData{:});
        currData = cellfun(@str2num, currData(:,1:9));

        % sort by first row and write to global storage
        currSortedData = sortrows(currData,1);
        sortedData{i} = currSortedData;
    end

    % save to mat file
    outName = ['tabTraj_nvt_pt' num2str(ff) '.mat'];
    save(outName,'sortedData','volData','-v7.3');
end
