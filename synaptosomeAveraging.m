
clear all
close all
clc

% Filename pattern
pattern = 'phys_BC';

% Get path to directory with images or ROIs
directory = '/Volumes/Ezra1/PhD in Chemistry/analysis_20190108/4Cb_results/Results_combined/Repeat_Oct-Nov/ripley/phys/images';

% Create new output folder
outputpath = fullfile(directory,'Particle_averaging');
mkdir(outputpath);

% Get list of files in directory with pattern in filename
if ispc
    filelist = dir([directory ['\*' pattern '*']]);
elseif ismac
    filelist = dir([directory ['/*' pattern '*']]);
end

% Loop over files
count = 0;
for i=1:length(filelist)
    
    % Get current filename from filelist
    filename = filelist(i).name;
    
    % Skip file if name starts with '.'
    if filename(1) == '.'
        continue;
    else
        count = count + 1;
    end
    
    % Get full path to current file
    filepath = fullfile(directory,filename);

    % Read image
    img = imread(filepath);
    
    % Sum with previous image
    if count > 1
        img_stack = cat(3,img_stack,img);
    else
        img_stack = img;
    end
    
end

% Get average image
img_avg = mean(img_stack,3);
img_avg = uint8(img_avg);
img_name = strcat(char(pattern),'_average.png');
imwrite(img_avg, fullfile(outputpath,img_name));

% Get median image
img_med = median(double(img_stack),3);
img_med = uint8(img_med);
img_name = strcat(char(pattern),'_median.png');
imwrite(img_med, fullfile(outputpath,img_name));

% Get standard deviation image
img_std = std(double(img_stack),0,3);
img_std = uint8(img_std);
img_name = strcat(char(pattern),'_stdev.png');
imwrite(img_std, fullfile(outputpath,img_name));

% Display
figure;
subplot(131)
imshow(img_avg);
title('Average')
subplot(132)
imshow(img_med);
title('Median')
subplot(133)
imshow(img_std);
title('Standard deviation')
