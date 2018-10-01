function [ PathName, FileNames_ch ] = GetImages( DefaultPath, channel_token )
% This function fetches the names of the red, green, and blue channel files for
% the .png visualizations of a reconstructed STORM image. 

% Author: Pedro Vallejo Ramirez ppv23@cam.ac.uk
% Laser Analytics Group
% Updated 28/09/2018


PathName = uigetdir(DefaultPath , 'Choose directory containing the png STORM images...');
FileList_ch = dir([PathName,filesep,'*',channel_token,'*','.png']); % List all the red channel files

N_files = length(FileList_ch); %number of files
FileNames_ch = cell(N_files,1);
    
        for i = 1:N_files
            FileNames_ch{i} = [PathName, filesep, FileList_ch(i).name];
        end
     FileNames_ch = natsortfiles(FileNames_ch);
end  