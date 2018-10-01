% This script imports individual synaptosome images and displaying them in
% a montage with their labels. The script looks for the cropped 8-bit
% visualizations of individual synaptosomes output by compareConditions.m
% with images for three different fluorescent channels, and displays them
% in a montage with an overlay of the three channels and the label for the
% synaptosome ID. 

% Ancillary functions:

% GetImages.m
% natsortfiles.m
% imdisp.m (Copyright (c) 2010, Oliver Woodford. All rights reserved.)

% Author: Pedro Vallejo Ramirez ppv23@cam.ac.uk
% Laser Analytics Group
% Updated 28/09/2018

%% Load raw synaptosome reconstructions 
clear all
close all
Path = 'E:\Experiments\synaptosomes\analysis_2ndRound\Results_combined\ripley\phys\images';

% select area token and label token for your data files

RC_token = '_phys_RC_synaptosome_';
GC_token = '_phys_GC_synaptosome_';
BC_token = '_phys_BC_synaptosome_';

[PathNameRC, FileNames_RC] = GetImages(Path, RC_token );
[PathNameGC, FileNames_GC] = GetImages(Path, GC_token );
[PathNameBC, FileNames_BC] = GetImages(Path, BC_token );

N_files = size(FileNames_RC,1);
disp('------------------------------------------------');
disp('------------------------------------------------');

disp('Opening files in following folder:');
disp(PathNameRC);
disp(['Number of files found: ', num2str(N_files)]);

%% display images 
imred = cell(N_files,1);
imgreen = cell(N_files,1);
imblue = cell(N_files,1);
overlay = cell(N_files,1);
red = cell(N_files,1);
green = cell(N_files,1);
blue = cell(N_files,1);
imlabel = cell(N_files,1);
im_names = cell(N_files,1);


for i = 1: N_files
    % read in images and transform into RGB to get a colormap for each
    imred{i,1}      = imread(FileNames_RC{i});
    imgreen{i,1}    = imread(FileNames_GC{i});
    imblue{i,1}     = imread(FileNames_BC{i});
    im_size         = size(imred{1,1});
    overlay{i,1}    = cat(3,imred{i,1},imgreen{i,1},imblue{i,1});    
    red{i,1}        = cat(3,imred{i,1},zeros(size(imred{i,1})),zeros(size(imred{i,1})));
    green{i,1}      = cat(3,zeros(size(imred{i,1})),imgreen{i,1},zeros(size(imred{i,1})));
    blue{i,1}       = cat(3,zeros(size(imred{i,1})),zeros(size(imred{i,1})),imblue{i,1});
    
    % write out labels for each image
    imlabel{i,1}    = ones(im_size(1),im_size(1));
    im_size2        = size(imlabel{i,1});
    [filepath,name,ext] = fileparts(FileNames_RC{i});
    string_length       = strlength(name);
    string              = strsplit(name,'_');
    im_names{i,1}       = [string{1} ' ' string{2} ' ' string{5}];
    position            = [im_size2(2)/2 im_size2(1)/2];
    imlabel{i,1}        = insertText(imlabel{i,1},position,im_names{i,1},'AnchorPoint','Center','TextColor','black','FontSize',16,'BoxColor','white');
    imlabel{i,1}        = rot90(imlabel{i,1});
    
end 

% Display image montage
figure;imdisp([red green blue overlay imlabel],'Size',[5 N_files],'Border',0.02)



