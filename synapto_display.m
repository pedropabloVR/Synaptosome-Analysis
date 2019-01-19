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
close all


Path        = 'E:\Experiments\synaptosomes\analysis_20190107\37C_results_thresh30\Results_combined\Repeat_B\ripley\phys\images';
output_dir  = 'E:\Experiments\synaptosomes\analysis_20190107\37C_results_thresh30\Results_combined';

% load .csv file with trimmed pooled results after manual removal of "ugly"
% regions

csv_path = "E:\Experiments\synaptosomes\analysis_20190107\4C_results\Results_combined\results_combined_pooled_trimmed.csv";
results_pooled = resultscombinedpooledtrimmed;

% select area token and label token for your data files
condition = 'phys';
repeat    = 'B';
im_cutoff = 20;
show      = 0;
RC_token  = ['_' condition '_RC_synaptosome_'];
GC_token  = ['_' condition '_GC_synaptosome_'];
BC_token  = ['_' condition '_BC_synaptosome_'];

[PathNameRC, FileNames_RC] = GetImages(Path, RC_token );
[PathNameGC, FileNames_GC] = GetImages(Path, GC_token );
[PathNameBC, FileNames_BC] = GetImages(Path, BC_token );

N_files = size(FileNames_RC,1);
disp('------------------------------------------------');
disp('------------------------------------------------');

disp('Opening files in following folder:');
disp(PathNameRC);
disp(['Number of files found: ', num2str(N_files)]);

% choose names from relevant repeat and condition
sampleIDs_repeat = results_pooled(strcmp(results_pooled.Repeat,repeat),:);
sampleIDs_condition = sampleIDs_repeat(strcmp(sampleIDs_repeat.condition,condition),:);
num_trimmed_files = size(sampleIDs_condition);
sampleIDs_names = cell(num_trimmed_files(1),1);
% concatenate names of sample ID and synaptosome ID
for i = 1:num_trimmed_files(1)
    temp = strsplit(sampleIDs_condition.sampleID(i),'_');
    samp_ID = strcat(temp(1)," ",temp(2));
    sampleIDs_names{i,1} = strcat(samp_ID," " ,num2str(sampleIDs_condition.synaptosomeID(i)));

end 

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
trimmed_names_indices = false(N_files,1);


path_output = fullfile(output_dir, filesep,['Synapto_montage_' condition repeat ]);
if exist(path_output, 'dir')
    opts.Interpreter = 'tex';
    opts.Default = 'Continue';
    quest = '\fontsize{12}An output folder ''Results'' already exists. If you continue, data in this folder might be overwritten.';
    answer = questdlg(quest,'Message','Cancel','Continue',opts);
    if strcmp(answer,'Continue')
        mkdir(path_output);
    else
        return
    end
else
    mkdir(path_output);
end


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
    
    % now choose only images with labels which correspond to labels on the
    % trimmed spreadsheet
    for lp = 1:num_trimmed_files(1)
        same = strcmp(im_names{i,1},sampleIDs_names{lp});
        if same 
         trimmed_names_indices(i) = true;  
        end    
    end 
    
end 

% Keep only images which correspond to selected ones on the trimmed
% spreadsheet
red     = red(trimmed_names_indices,:);
green   = green(trimmed_names_indices,:);
blue    = blue(trimmed_names_indices,:);
overlay = overlay(trimmed_names_indices,:);
imlabel = imlabel(trimmed_names_indices,:);
N_files = size(red,1);
% cycle through the images, 20 at a time, display them, and then continue
% with the next 20 until there are none left

% divide arrays into sub-arrays with 20 images each
remainder = mod(N_files,im_cutoff);
if remainder >0
    num_sub_arrays = floor(N_files/im_cutoff)+1;
else
    num_sub_arrays = floor(N_files/im_cutoff);
end 

count = 1;
if num_sub_arrays == 1
    remainder = N_files;
end 
for j = 1:num_sub_arrays
    % only on the last count of the loop, if the mod is more than 0, fix
    if j == num_sub_arrays % on the last run of the loop, run until the end of the number of files 
        red_j(1:remainder,1)   = red(count:N_files,1);
        green_j(1:remainder,1) = green(count:N_files,1);
        blue_j(1:remainder,1)  = blue(count:N_files,1);
        overlay_j(1:remainder,1)  = overlay(count:N_files,1);
        imlabel_j(1:remainder,1)  = imlabel(count:N_files,1);
        h=figure;imdisp([red_j green_j blue_j overlay_j imlabel_j],'Size',[5 remainder],'Border',0.02)

    else
        red_j(1:im_cutoff,1)   = red(count:j*im_cutoff,1);
        green_j(1:im_cutoff,1) = green(count:j*im_cutoff,1);
        blue_j(1:im_cutoff,1)  = blue(count:j*im_cutoff,1);
        overlay_j(1:im_cutoff,1)  = overlay(count:j*im_cutoff,1);
        imlabel_j(1:im_cutoff,1)  = imlabel(count:j*im_cutoff,1);
        count = count+im_cutoff;
        h=figure;imdisp([red_j green_j blue_j overlay_j imlabel_j],'Size',[5 im_cutoff],'Border',0.02)
    end 
    
    % export figure to a pdf for display
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    filename = [condition '_synapto_montage_' num2str(j)];
    print(h,[path_output filesep filename],'-dpdf','-r300')
    clear red_j green_j blue_j overlay_j imlabel_j 
end 

% Display image montage
%figure;imdisp([red green blue overlay imlabel],'Size',[5 N_files],'Border',0.02)
% 
% for j = 1:num_sub_arrays
%     % only on the last count of the loop, if the mod is more than 0, fix
%     if j == num_sub_arrays % on the last run of the loop, run until the end of the number of files
%         red_j(count:N_files,1) = red(count:N_files,1);
%         green_j(count:N_files,1) = green(count:N_files,1);
%         blue_j(count:N_files,1) = blue(count:N_files,1);
%     else
%         red_j(count:j*im_cutoff,1) = red(count:j*im_cutoff,1);
%         green_j(count:j*im_cutoff,1) = green(count:j*im_cutoff,1);
%         blue_j(count:j*im_cutoff,1) = blue(count:j*im_cutoff,1);
%         count = count+im_cutoff;
%     end 
%     
%     figure;imdisp([red_j green_j blue_j ],'Size',[3 im_cutoff],'Border',0.02)
% end 
% 


