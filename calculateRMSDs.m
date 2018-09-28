% This function reads in all localisation files (regions that contain only
% localisations of one cluster, e.g. regions cropped by synaptoCROP) in a
% folder, extracts the x and y coordinates and calculated the RMSD. The
% function returns an array containing the rmsd calculated for the
% individual localisation files in the folder, as well as the mean and
% standard deviation.
% 
% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated on 19 Sept 2018


function [rmsd, mean_rmsd, std_rmsd] = calculateRMSDs(path_folder,format)

% Get list of files in folder
if strcmp(format,'thunderstorm')
    filelist = struct2cell(dir([path_folder,filesep,'*.csv']));
elseif strcmp(format,'rapidstorm')
    filelist = struct2cell(dir([path_folder,filesep,'*.txt']));
else
    disp(['The format ' format ' doesn''t exist. Please pick a valid format in calculateRMSDs.']);
end
filelist = filelist(1,:);
number_of_files = length(filelist);

disp(['Folder: ' path_folder]);
disp(['Number of files found: ', num2str(number_of_files)]);

% Read in the files and calculate RMSD
rmsd = zeros(1,number_of_files);
for i = 1:number_of_files
    % Get filename and path
    filename = filelist{i};
    path = fullfile(path_folder,filename);
    % Read localisation file
    locs = readLocFile(path,format);
    rmsd(i) = getRMSD(locs.x,locs.y);
end
mean_rmsd = mean(rmsd);
std_rmsd = std(rmsd);
disp(['RMSD: ' num2str(mean_rmsd) ' ± ' num2str(std_rmsd) ' nm  (mean ± sd)']);
disp(' ')

end