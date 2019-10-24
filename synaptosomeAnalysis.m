%% Synapto-Analysis

% This script reads in a set of reconstructed single molecule localisation
% data sets from three different color channels and it automatically
% detects individual synaptosomes from each reconstructed field of view. 
% The script performs simple filtering on the data sets based on frame
% number, sigma value of the Gaussian fit, localisation precision, and
% intensity. A density filter is implemented to remove background
% localisations and only keep clustered detections. 
%
% The script outputs:
% - Reconstructed images for the three different color channels 
% - Overlay image of the three channels 
% - Binary image with potential detected synaptosomes in the red channel 
% - Text file with a summary of the thresholds established, and potential
%   number of synaptosomes detected. 
% - CSV file with centroids and weighted overlap calculations for each
%   potentially detected synaptosome. 

% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated: 13 September 2019


clear all
close all
clc

tic

%% File parameters

%directory = 'F:\synaptosomes\2018_10_10_Pedro_5thRound_EGTAK\output_reconstructions\Registered_data';
%directory = 'E:\Experiments\synaptosomes\Datasets_synaptosomes_20181206_4C_37C\37Cb\Data\thunderSTORM_phys\reconstructions\Registered_data';
%directory ='F:\synaptosomes\2018_10_10_Pedro_5thRound_PHYS_ISR\output_reconstructions\Registered_data'; % last setting used before 12/09/2019
directory = 'D:\Experiments\synaptosomes\Datasets_synaptosomes_20181206_4C_37C\4Cb\Data\thunderSTORM_phys\reconstructions\Registered_data';
% path to folder where outputfolder will be created (if doesn't already exist)
output_dir = fullfile('D:\Experiments\synaptosomes\test for supplement figures 20191023',filesep);

% remember to change the format for filtering to rapidSTORM when using the
% 37C data sets. 
repeat           = '4Cb';
condition        = 'phys';

channel_token_RC = '_647';
channel_token_GC = '_561_reg';
channel_token_BC = '_488_reg';
rotate_channel2  = 0;

% Single molecule reconstruction settings

pixelsize            = 117; % pixelsize in nm
magnification        = 10; % value of 10 gives 11.7 nm pixels in reconstruction (if pixelsize camera is 117 nm)
show                 = 1; % 1 to show extra intermediate results
show_d               = 1;
format               = 'thunderstorm'; % reconstruction software used (only thunderstorm)
format_for_filtering = 'thunderstorm'; % SWITCH TO RAPIDSTORM FOR THE 37C DATA

% Filtering parameters
% set these flags to 1 to activate filters
density_filter  = 1;
filter          = 1;    % to filter the localisation files before further processing based on the following parameters

% for the density filters
max_radius_RC  = 1500; % maximum radius of clusters in red channel
min_nr_locs_RC = 300;  % minimum nr of locs within min_radius_RC from each localisation in red channel for density filtering
max_radius_GC  = 0;    % maximum radius of clusters in green channel
min_nr_locs_GC = 0;    % minimum nr of locs within min_radius_GC from each localisation in green channel
max_radius_BC  = 0;    % maximum radius of clusters in blue channel
min_nr_locs_BC = 0;    % minimum nr of locs within min_radius_BC from each localisation in blue channel

P_RC           = 200; % blobs with fewer pixels will be removed from mask of red channel
P_GC           = 0;   % blobs with fewer pixels will be removed from mask of green channel
P_BC           = 0;   % blobs with fewer pixels will be removed from mask of blue channel

level_RC       = 0; % threshold (0-1) for mask red   channel (0 for Otsu's threshold)
level_GC       = 0; % threshold (0-1) for mask green channel (0 for Otsu's threshold)
level_BC       = 0; % threshold (0-1) for mask blue  channel (0 for Otsu's threshold)

sigma_kernel   = 15; % sigma for generating an image (if format = rapidstorm)
remove_border  = 0; % to remove the edges of images if there are edge artefacts
border_width   = 700; % width of border removed (in image coordinates)

% for the localisation filters
min_frames      = 500;  % throw away initial frames to avoid really high density data
min_sigma       = 40;   % maximum sigma value of the Gaussian fit
max_sigma       = 300;  % maximum sigma value of the Gaussian fit
min_intensity   = 500;  % min number of photons expected for an event
max_uncertainty = 40;   % max uncertainty expected for single molecules in focus
min_uncertainty = 5;    % max uncertainty expected for single molecules in focus 

essence = 0; % 1 to only write away files essential for further analysis
warning('off','images:initSize:adjustingMag');
warning('off','MATLAB:MKDIR:DirectoryExists');


%% Creating output folder and parameter file

% Create new output folder
path_output = fullfile(output_dir,['Results_' repeat ' ' condition]);
if exist(path_output, 'dir')
    opts.Interpreter = 'tex';
    opts.Default = 'Continue';
    quest = '\fontsize{12}An output folder ''Results'' already exists. If you continue, data in this folder might be overwritten.';
    answer = questdlg(quest,'Message','Cancel','Continue',opts);
    if strcmp(answer,'Continue')
        mkdir(path_output);
        mkdir(fullfile(path_output,'regions'));
    else
        return
    end
else
    mkdir(path_output);
    mkdir(fullfile(path_output,'regions'));
end


% Write away workspace variables in a parameter file
save(fullfile(path_output,'parameters.mat'));
 max_radius_RC = max_radius_RC/magnification;
 max_radius_GC = max_radius_GC/magnification;
 max_radius_BC = max_radius_BC/magnification;

% Get list of locfiles of red channel
if ismac
    filelist = dir([directory strcat('/*',channel_token_RC,'.csv')]);
elseif ispc 
    filelist = dir([directory strcat('\*',channel_token_RC,'.csv')]);
end
    
for i = 1:size(filelist,1)
    
    % Get filenames
    filename_RC = filelist(i).name;
    area_token = strsplit(filename_RC,'_'); area_token = area_token{1};
    filename_GC = strcat(area_token,channel_token_GC,'.csv');
    filename_BC = strcat(area_token,channel_token_BC,'.csv');
    
    % Write some info to command window
    disp('##### Synaptosome analysis #####')
    disp(' ')
    disp(['Folder:             ' directory])
    disp(['File red channel:   ' filename_RC])
    disp(['File green channel: ' filename_GC])
    disp(['File blue channel:  ' filename_BC])
    
    % Open a summary file for writing
    summary_file = fopen(fullfile(path_output,strcat(area_token,'_summary.txt')),'wt');
    
    % Read in files
    locs_RC = readLocFile(fullfile(directory,filename_RC),format);
    locs_GC = readLocFile(fullfile(directory,filename_GC),format);
    locs_BC = readLocFile(fullfile(directory,filename_BC),format);
    
    % Record number of localisations before filtering
    num_locs_RC_prefilter = size(locs_RC,1);
    num_locs_GC_prefilter = size(locs_GC,1);
    num_locs_BC_prefilter = size(locs_BC,1);
    
%     % Write away raw data from localisations
%     
%     writeLocFile(locs_RC,fullfile(path_output,strcat(area_token,'_locs_RC_raw.csv')),format)
%     writeLocFile(locs_GC,fullfile(path_output,strcat(area_token,'_locs_GC_raw.csv')),format)
%     writeLocFile(locs_BC,fullfile(path_output,strcat(area_token,'_locs_BC_raw.csv')),format)
%       
    
    %% Estimate Fov
    % Assumes there will be at least one localisation close to the edges
    % of the Fov in one of the three channels
    
    Fov = estimateFov([locs_RC.x; locs_GC.x; locs_GC.x; locs_RC.y; locs_BC.y; locs_BC.y], pixelsize);
    disp(['Estimated Fov:      ' num2str(Fov)])
    disp(' ')
    
    
    %% Optional removal of edge artefacts
    % Big high intensity blobs around the edges of images can throw off the
    % synaptosome detector, so a border with some specified width can be
    % cropped out before doing any further processing.
    
    % For removing edge artefacts
    if remove_border
        %   NOTE: NEED TO INCLUDE CODE TO REMOVE Y POSITIONS FROM THE
        %   TRIMMED X DATA, AND TO REMOVE X POSITIONS FROM THE TRIMMED Y
        %   DATA. 
        locs_RC_filtered = locs_RC(locs_RC.x > border_width,:);
        locs_GC_filtered = locs_GC(locs_GC.x > border_width,:);
        locs_BC_filtered = locs_BC(locs_BC.x > border_width,:);
        
        locs_RC_filtered = locs_RC_filtered(locs_RC_filtered.x < Fov*pixelsize - border_width,:);
        locs_GC_filtered = locs_GC_filtered(locs_GC_filtered.x < Fov*pixelsize - border_width,:);
        locs_BC_filtered = locs_BC_filtered(locs_BC_filtered.x < Fov*pixelsize - border_width,:);
        
        locs_RC_filtered = locs_RC_filtered(locs_RC_filtered.y > border_width,:);
        locs_GC_filtered = locs_GC_filtered(locs_GC_filtered.y > border_width,:);
        locs_BC_filtered = locs_BC_filtered(locs_BC_filtered.y > border_width,:);
        
        locs_RC_filtered = locs_RC_filtered(locs_RC_filtered.y < Fov*pixelsize - border_width,:);
        locs_GC_filtered = locs_GC_filtered(locs_GC_filtered.y < Fov*pixelsize - border_width,:);
        locs_BC_filtered = locs_BC_filtered(locs_BC_filtered.y < Fov*pixelsize - border_width,:);
    else
        locs_RC_filtered = locs_RC;
        locs_GC_filtered = locs_GC;
        locs_BC_filtered = locs_BC;
    end
    
    
    %% Optional filtering
    % The registered localisations filtered to keep only localisations
    % detected after 500 frames (because these early frames are usually not
    % sparse enough to assure good reconstruction), keep only localisations
    % that have a sigma between 40 nm and 400 nm (to avoid grid artefacts
    % and exclude localizations that are out of focus) and uncertainty
    % smaller than 40 nm and larger than 5 nm (the expected resolution of the reconstructed
    % dSTORM images). For the red channel the intensity threshold is set at
    % 1000 photons, and for the green and blue channels at 500 photons. 
    
    if filter
        % Filtering red channel ---------------------------------------------------
        locs_RC_filtered = locs_RC_filtered(locs_RC_filtered.frame       > min_frames ,:);
        num_locs_RC_framefilter = size(locs_RC_filtered,1);
        if ~strcmp(format_for_filtering,'rapidstorm')
            locs_RC_filtered = locs_RC_filtered(locs_RC_filtered.sigma   > min_sigma  ,:);
            locs_RC_filtered = locs_RC_filtered(locs_RC_filtered.sigma   < max_sigma ,:);
        end
        locs_RC_filtered = locs_RC_filtered(locs_RC_filtered.intensity   > min_intensity,:);
        locs_RC_filtered = locs_RC_filtered(locs_RC_filtered.uncertainty < max_uncertainty  ,:);
        locs_RC_filtered = locs_RC_filtered(locs_RC_filtered.uncertainty > min_uncertainty  ,:);
        
        % Filtering green channel -------------------------------------------------
        locs_GC_filtered = locs_GC_filtered(locs_GC_filtered.frame       > min_frames,:);
        num_locs_GC_framefilter = size(locs_GC_filtered,1);
        if ~strcmp(format_for_filtering,'rapidstorm')
            locs_GC_filtered = locs_GC_filtered(locs_GC_filtered.sigma   > min_sigma  ,:);
            locs_GC_filtered = locs_GC_filtered(locs_GC_filtered.sigma   < max_sigma,:);
        end
        
        locs_GC_filtered = locs_GC_filtered(locs_GC_filtered.intensity   > min_intensity,:);
        locs_GC_filtered = locs_GC_filtered(locs_GC_filtered.uncertainty < max_uncertainty ,:);
        locs_GC_filtered = locs_GC_filtered(locs_GC_filtered.uncertainty > min_uncertainty ,:);
        
        % Filtering blue channel --------------------------------------------------
        locs_BC_filtered = locs_BC_filtered(locs_BC_filtered.frame       > min_frames,:);
        num_locs_BC_framefilter = size(locs_BC_filtered,1);
        if ~strcmp(format_for_filtering,'rapidstorm')
            locs_BC_filtered = locs_BC_filtered(locs_BC_filtered.sigma   > min_sigma  ,:);
            locs_BC_filtered = locs_BC_filtered(locs_BC_filtered.sigma   < max_sigma,:);
        end      
        locs_BC_filtered = locs_BC_filtered(locs_BC_filtered.intensity   > min_intensity,:);
        locs_BC_filtered = locs_BC_filtered(locs_BC_filtered.uncertainty < max_uncertainty ,:);
        locs_BC_filtered = locs_BC_filtered(locs_BC_filtered.uncertainty > min_uncertainty ,:);
    end
    
    % Get number of localisations after filter
    num_locs_RC_postfilter = size(locs_RC_filtered,1);
    num_locs_GC_postfilter = size(locs_GC_filtered,1);
    num_locs_BC_postfilter = size(locs_BC_filtered,1);
    
    % calculate overall rejection rate (percentage) after localisation filters
    % currently the base number of localisations is calculated after the
    % first 500 frames have been removed (the pre-filtered number of
    % localisations). 
    if filter
    overall_rejection_RC = ((num_locs_RC_framefilter-num_locs_RC_postfilter)/num_locs_RC_framefilter)*100;
    overall_rejection_GC = ((num_locs_GC_framefilter-num_locs_GC_postfilter)/num_locs_GC_framefilter)*100;
    overall_rejection_BC = ((num_locs_BC_framefilter-num_locs_BC_postfilter)/num_locs_BC_framefilter)*100;
    
    fprintf(summary_file,'Overall rejection rate after localisation filtering \n');
    fprintf(summary_file,'Red   channel: %2.f %%\n',overall_rejection_RC);
    fprintf(summary_file,'Green channel: %2.f %%\n',overall_rejection_GC);
    fprintf(summary_file,'Green channel: %2.f %%\n',overall_rejection_BC);
    end 
    %% Write away filtered localisations
    
    writeLocFile(locs_RC_filtered,fullfile(path_output,strcat(area_token,'_locs_RC_filtered.csv')),format)
    writeLocFile(locs_GC_filtered,fullfile(path_output,strcat(area_token,'_locs_GC_filtered.csv')),format)
    writeLocFile(locs_BC_filtered,fullfile(path_output,strcat(area_token,'_locs_BC_filtered.csv')),format)
       
    %% Generate SMLM reconstructions (before density filtering)
    
    % Red channel
    X_RC = locs_RC_filtered.x;
    Y_RC = locs_RC_filtered.y;
    if strcmp(format_for_filtering,'thunderstorm')
        if (locs_RC_filtered.sigma == 0)
            sigma_RC = sigma_kernel;
        else 
            sigma_RC = median(locs_RC_filtered.sigma)/magnification;
        end 
    elseif strcmp(format_for_filtering,'rapidstorm')
        sigma_RC = sigma_kernel;
    end
    %sigma_RC = sigma_kernel;
    intensities = locs_RC_filtered.intensity;
    img_RC = generateImage(X_RC,Y_RC,sigma_RC,intensities,Fov,pixelsize,magnification);
    filename = char(strcat(area_token,'_reconstruction_RC.tif'));
    if ~essence; imwrite(flip(flip(img_RC,2),1), fullfile(path_output, filename)); end
    if show; figure; imshow(flip(flip(img_RC,2),1)); title('Reconstruction red channel'); end
    
    % Green channel
    X_GC = locs_GC_filtered.x;
    Y_GC = locs_GC_filtered.y;
    if strcmp(format,'thunderstorm')
         if (locs_GC_filtered.sigma == 0)
            sigma_GC = sigma_kernel;
         else
            sigma_GC = median(locs_GC_filtered.sigma)/magnification;
        end 
    elseif strcmp(format,'rapidstorm')
        sigma_GC = sigma_kernel;
    end
    %sigma_GC = sigma_kernel;
    intensities = locs_GC_filtered.intensity;
    img_GC = generateImage(X_GC,Y_GC,sigma_GC,intensities,Fov,pixelsize,magnification);
    filename = char(strcat(area_token,'_reconstruction_GC.tif'));
    if ~essence; imwrite(flip(flip(img_GC,2),1), fullfile(path_output, filename)); end
    if show; figure; imshow(flip(flip(img_GC,2),1)); title('Reconstruction green channel'); end
    
    % Blue channel
    X_BC = locs_BC_filtered.x;
    Y_BC = locs_BC_filtered.y;
    if strcmp(format,'thunderstorm')
         if (locs_BC_filtered.sigma == 0)
            sigma_BC = sigma_kernel;
         else
            sigma_BC = median(locs_BC_filtered.sigma)/magnification;
        end 
    elseif strcmp(format,'rapidstorm')
        sigma_BC = sigma_kernel;
    end
    %sigma_BC = sigma_kernel;
    intensities = locs_BC_filtered.intensity;
    img_BC = generateImage(X_BC,Y_BC,sigma_BC,intensities,Fov,pixelsize,magnification);
    filename = char(strcat(area_token,'_reconstruction_BC.tif'));
    if ~essence; imwrite(flip(flip(img_BC,2),1), fullfile(path_output, filename)); end
    if show; figure; imshow(flip(flip(img_BC,2),1)); title('Reconstruction blue channel'); end
    
    % Get two-colour reconstruction red-green (only for visualization)
    merged_reconstruction_RG = uint8(zeros(size(img_RC,1),size(img_RC,2),3));
    merged_reconstruction_RG(:,:,1) = uint8(img_RC);
    merged_reconstruction_RG(:,:,2) = uint8(img_GC);
    filename = char(strcat(area_token,'_reconstruction_RG.tif'));
    if ~essence; imwrite(flip(flip(merged_reconstruction_RG,2),1), fullfile(path_output, filename)); end
    if show; figure; imshow(flip(flip(merged_reconstruction_RG,2),1)); title('Reconstruction RG'); end
    
    % Get two-colour reconstruction red-blue (only for visualization)
    merged_reconstruction_RB = uint8(zeros(size(img_RC,1),size(img_RC,2),3));
    merged_reconstruction_RB(:,:,1) = uint8(img_RC);
    merged_reconstruction_RB(:,:,3) = uint8(img_BC);
    filename = char(strcat(area_token,'_reconstruction_RB.tif'));
    if ~essence; imwrite(flip(flip(merged_reconstruction_RB,2),1), fullfile(path_output, filename)); end
    if show; figure; imshow(flip(flip(merged_reconstruction_RB,2),1)); title('Reconstruction RB'); end
    
    % Get three-colour reconstruction (only for visualization)
    merged_reconstruction_RGB = uint8(zeros(size(img_RC,1),size(img_RC,2),3));
    merged_reconstruction_RGB(:,:,1) = uint8(img_RC);
    merged_reconstruction_RGB(:,:,2) = uint8(img_GC);
    merged_reconstruction_RGB(:,:,3) = uint8(img_BC);
    filename = char(strcat(area_token,'_reconstruction_RGB.tif'));
    imwrite(flip(flip(merged_reconstruction_RGB,2),1), fullfile(path_output, filename));
    if show; figure; imshow(flip(flip(merged_reconstruction_RGB,2),1)); title('Reconstruction RGB'); end
    
    
    %% Density filtering
    % To improve detection of synaptosomes in the red channel, some density
    % filtering is applied (the images are cluttered). Only localisations
    % that have a minimum number of neighbouring localisations within some
    % specified search radius are kept. Considering the radius of a
    % synaptosome is about 500 nm, this value is chosen for the search
    % radius. By looking at a number of images and testing out different
    % parameters, the minimum number of localisations within this radius
    % was chosen to be 500. The variables containing the unfiltered
    % localisations are cleared after filtering to save memory.
    
    fprintf(summary_file,'Rejection Rate after density filtering \n');
    if ~density_filter
            fprintf(summary_file,'Density filter off \n');
    end 
    if density_filter
        % Density filtering red channel -------------------------------------------
        if max_radius_RC ~= 0 && min_nr_locs_RC ~= 0
            % Perform density filtering
            [locs_RC_density_filtered,indeces_RC,~] = ...
                nearestNeighbourDensityFilter([locs_RC_filtered.x locs_RC_filtered.y], ...
                max_radius_RC,min_nr_locs_RC,show_d); disp(' ');
            if isempty(locs_RC_density_filtered)
                disp('All localizations in the red channel were filtered out during density filtering!');
                continue
            else
                nr_locs_before = size(locs_RC_filtered.x,1);
                nr_locs_after  = size(locs_RC_density_filtered,1);
                fprintf(summary_file,'Red   channel: %1.f %%\n', ((nr_locs_before-nr_locs_after)/nr_locs_before)*100);
            end
            % Generate image from filtered localizations
            X_RC = locs_RC_density_filtered(:,1);
            Y_RC = locs_RC_density_filtered(:,2);
            if strcmp(format_for_filtering,'thunderstorm')
                sigma_RC = median(locs_RC_filtered.sigma)/magnification;
            elseif strcmp(format_for_filtering,'rapidstorm')
                sigma_RC = sigma_kernel;
            end
            %sigma_RC = sigma_kernel;
            intensities = array2table([indeces_RC locs_RC_filtered.intensity],'VariableNames',{'id','intensity'});
            intensities_filtered = intensities(intensities.id == 1,:);
            intensities_filtered = intensities_filtered.intensity;

            img_RC = generateImage(X_RC,Y_RC,sigma_RC,intensities_filtered,Fov,pixelsize,magnification);
            filename = char(strcat(area_token,'_reconstruction_RC_density_filtered.tif'));
            if ~essence; imwrite(flip(flip(img_RC,2),1), fullfile(path_output, filename)); end
            if show; figure; imshow(flip(flip(img_RC,2),1)); title('Reconstruction red channel after density filtering'); end
        else
            locs_RC_density_filtered = [locs_RC_filtered.x locs_RC_filtered.y];
            fprintf(summary_file,'Red   channel: 0 %%\n');
        end
        clear var locs_RC % Clear variable containing unfiltered data to save memory

        % Density filtering green channel -----------------------------------------
        if max_radius_GC ~= 0 && min_nr_locs_GC ~= 0
            [locs_GC_density_filtered, indeces_GC, numNeighbours_GC] = ...
                nearestNeighbourDensityFilter([locs_GC_filtered.x locs_GC_filtered.y], ...
                max_radius_GC,min_nr_locs_GC,show); disp(' ');
            if isempty(locs_GC_density_filtered) % If no localizations remain after density filtering
                disp('All localizations in the green channel were filtered out during density filtering!');
                return
            else
                nr_locs_before = size(locs_GC_filtered.x,1);
                nr_locs_after  = size(locs_GC_density_filtered,1);
                fprintf(summary_file,'Green channel: %1.f %%\n', ((nr_locs_before-nr_locs_after)/nr_locs_before)*100);
            end
            % Generate image from filtered localisations
            X_GC = locs_GC_density_filtered(:,1);
            Y_GC = locs_GC_density_filtered(:,2);
            if strcmp(format_for_filtering,'thunderstorm')
                sigma_GC = median(locs_GC_filtered.sigma)/magnification;
            elseif strcmp(format_for_filtering,'rapidstorm')
                sigma_GC = sigma_kernel;
            end
            %sigma_GC = sigma_kernel;
            intensities = array2table([indeces_GC locs_GC_filtered.intensity],'VariableNames',{'id','intensity'});
            intensities_filtered = intensities(intensities.id == 1,:);
            intensities_filtered = intensities_filtered.intensity;

            img_GC = generateImage(X_GC,Y_GC,sigma_GC,intensities_filtered,Fov,pixelsize,magnification);
            filename = char(strcat(area_token,'_reconstruction_GC_density_filtered.tif'));
            if ~essence; imwrite(flip(flip(img_GC,2),1), fullfile(path_output, filename)); end
            if show; figure; imshow(flip(flip(img_GC,2),1)); title('Reconstruction green channel after density filtering'); end
        else
            locs_GC_density_filtered = [locs_GC_filtered.x locs_GC_filtered.y];
            fprintf(summary_file,'Green channel: 0 %%\n');
        end
        clear var locs_GC % Clear variable containing unfiltered data to save memory

        % Density filtering blue channel ------------------------------------------
        if max_radius_BC ~= 0 && min_nr_locs_BC ~= 0
            [locs_BC_density_filtered, indeces_BC, numNeighbours_BC] = ...
                nearestNeighbourDensityFilter([locs_BC_filtered.x locs_BC_filtered.y], ...
                max_radius_BC,min_nr_locs_BC,show); disp(' ');
            if isempty(locs_BC_density_filtered) % If no localizations remain after density filtering
                disp('All localizations in the blue channel were filtered out during density filtering!');
                return
            else
                nr_locs_before = size(locs_BC_filtered.x,1);
                nr_locs_after  = size(locs_BC_density_filtered,1);
                fprintf(summary_file,'Blue  channel: %1.f %%\n\n', ((nr_locs_before-nr_locs_after)/nr_locs_before)*100);
            end
            % Generate image from filtered localizations
            X_BC = locs_BC_density_filtered(:,1);
            Y_BC = locs_BC_density_filtered(:,2);
            if strcmp(format_for_filtering,'thunderstorm')
                sigma_BC = median(locs_BC_filtered.sigma)/magnification;
            elseif strcmp(format_for_filtering,'rapidstorm')
                sigma_BC = sigma_kernel;
            end
            %sigma_BC = sigma_kernel;
            intensities = array2table([indeces_BC locs_BC_filtered.intensity],'VariableNames',{'id','intensity'});
            intensities_filtered = intensities(intensities.id == 1,:);
            intensities_filtered = intensities_filtered.intensity;

            img_BC = generateImage(X_BC,Y_BC,sigma_BC,intensities_filtered,Fov,pixelsize,magnification);
            filename = char(strcat(area_token,'_reconstruction_BC_density_filtered.tif'));
            if ~essence; imwrite(flip(flip(img_BC,2),1), fullfile(path_output, filename)); end
            if show; figure; imshow(flip(flip(img_BC,2),1)); title('Reconstruction blue channel after density filtering'); end
        else
            locs_BC_density_filtered = [locs_BC_filtered.x locs_BC_filtered.y];
            fprintf(summary_file,'Blue  channel: 0 %%\n\n');
        end
        clear var locs_BC % Clear variable containing unfiltered data to conserve memory
    end   
    
    if rotate_channel2
        img_GC = rot90(img_GC);
    end 
        
    %% Get binary mask for all three channels
    % First an image is generated from the filtered localisations. The
    % magnification should be chosen to match the resolution of the
    % reconstruction (in our case the pixelsize of the camera is 117 nm and
    % the mean resolution of the reconstructions 30 nm (as determined by
    % FRC), so a magnification = 10 is used to get effective pixelsize of
    % 11.7 nm which meets the Nyquist criterion). The generated image is
    % thresholded using Otsu's threshold to get a binary mask. In the red
    % channel, blobs in the mask that are smaller than some amount of
    % pixels P are filled (because they probably are not synaptosomes).
    
    fprintf(summary_file,'Threshold (Otsu''s method or user-specified):\n');
    
    % Red channel -------------------------------------------------------------
    % Binarize image to get a mask
    if level_RC == 0 % Get automatic threshold
        level_RC_i = graythresh(img_RC);
        disp(['Automatic threshold for red   channel: ' num2str(level_RC_i)]);
    else
        level_RC_i = level_RC;
        disp(['Specified threshold for red   channel: ' num2str(level_RC_i)]);
    end
    fprintf(summary_file,'Red   channel: %f\n',level_RC_i);
    mask_RC = im2bw(img_RC,level_RC_i);
    
    % Remove blobs that are too small
    mask_RC = bwareaopen(mask_RC,P_RC);
    if show
        figure;
        imshow(flip(flip(mask_RC,2),1));
        title(sprintf('Red channel after bwareaopen (P = %d)', P_RC));
    end
    if ~essence
        filename = strcat(area_token,'_detected_synaptosomes.tif');
        imwrite(flip(flip(mask_RC,2),1), fullfile(path_output,filename));
    end
    
    
    % Green channel -----------------------------------------------------------
    % Binarize image to get a mask
    if level_GC == 0 % Get automatic threshold
        level_GC_i = graythresh(img_GC);
        disp(['Automatic threshold for green channel: ' num2str(level_GC_i)]);
    else
        level_GC_i = level_GC;
        disp(['Automatic threshold for green channel: ' num2str(level_GC_i)]);
    end
    fprintf(summary_file,'Green channel: %f\n',level_GC_i);
    mask_GC = im2bw(img_GC,level_GC_i);
    
    % Remove blobs that are too small
    mask_GC = bwareaopen(mask_GC,P_GC);
    if show
        figure;
        imshow(flip(flip(mask_GC,2),1));
        title(sprintf('Green channel after bwareaopen (P = %d)', P_GC));
    end
    
    % Blue channel ------------------------------------------------------------
    % Binarize image to get a mask
    if level_BC == 0 % Get automatic threshold
        level_BC_i = graythresh(img_BC);
        disp(['Automatic threshold for blue  channel: ' num2str(level_BC_i)]); disp(' ');
    else
        level_BC_i = level_BC;
        disp(['Automatic threshold for blue  channel: ' num2str(level_BC_i)]); disp(' ');
    end
    fprintf(summary_file,'Blue  channel: %f\n',level_BC_i);
    mask_BC = im2bw(img_BC,level_BC_i);
    
    % Remove blobs that are too small
    mask_BC = bwareaopen(mask_BC,P_BC);
    if show
        figure;
        imshow(flip(flip(mask_BC,2),1));
        title(sprintf('Blue channel after bwareaopen (P = %d)', P_BC));
    end
    
    % Get RGB image mask ------------------------------------------------------
    merged_mask = zeros(size(mask_RC,1),size(mask_RC,2),3);
    merged_mask(:,:,1) = mask_RC;
    merged_mask(:,:,2) = mask_GC;
    merged_mask(:,:,3) = mask_BC;
    filename = char(strcat(area_token,'_masksRGB.tif'));
    imwrite(flip(flip(merged_mask,2),1), fullfile(path_output, filename));
    if show; figure; imshow(flip(flip(merged_mask,2),1)); title('Merged channels'); end
    
    %% Calculate overlap between channels for all potential synaptosomes
    
    % Get labeled mask of red channel
    mask_RC_labeled = bwlabel(mask_RC);
        
    % Get number of potential synaptosomes detected
    numberOfClusters = max(max(mask_RC_labeled));
    disp([num2str(numberOfClusters) ' potential synaptosomes detected.']); disp(' ');
    fprintf(summary_file,'%d potential synaptosomes detected.\n',numberOfClusters);
    
    if numberOfClusters > 0
    
        % Get centroids of the synaptosomes
        synaptosomeMeasurements = regionprops(mask_RC_labeled,mask_RC_labeled,'all');
        centroids = [synaptosomeMeasurements.Centroid];
        x_centroid = centroids(1:2:end-1);
        y_centroid = centroids(2:2:end);

        % Get area of the synaptosomes
        areas = [synaptosomeMeasurements.Area];
        % areas = areas*(magnification^2); % correct here for magnification factor, or do it later in analysis

        % Loop over detected synaptosomes and calculate their overlap with clusters
        % in the green and blue channel as a percentage, and a weighted overlap
        % using the image intensities
        overlap_perc_GC = zeros(1,numberOfClusters);
        overlap_perc_BC = zeros(1,numberOfClusters);

        weighted_overlap_GC = zeros(1,numberOfClusters);
        weighted_overlap_BC = zeros(1,numberOfClusters);
        weighted_overlap_GB = zeros(1,numberOfClusters);
        
        
        overlap_RonG = zeros(1,numberOfClusters);
        overlap_GonR = zeros(1,numberOfClusters);
        overlap_RonB = zeros(1,numberOfClusters);
        overlap_BonR = zeros(1,numberOfClusters);
        overlap_GonB = zeros(1,numberOfClusters);
        overlap_BonG = zeros(1,numberOfClusters);
        
        
       % initialize cells to store cropped images used to calculate
       % overlap. 
        red_crop     = cell(1,numberOfClusters);
        green_crop   = cell(1,numberOfClusters);
        blue_crop    = cell(1,numberOfClusters);
        
        for j = 1:numberOfClusters

            disp(['Synaptosome ' num2str(j) '/' num2str(numberOfClusters)]);

            % Get mask for current synaptosome only
            maskCurrentSynaptosome = mask_RC_labeled;
            maskCurrentSynaptosome(maskCurrentSynaptosome~=j) = 0;
            maskCurrentSynaptosome(maskCurrentSynaptosome==j) = 1;

            % make dilated mask of synaptosome location and obtain bounding box
            sr = strel('square',45);
            maskBoundingRegion = imdilate(maskCurrentSynaptosome,sr);
            boundingregionprops = regionprops(maskBoundingRegion);
            
            % use bounding box to extract 
            %mask_boundingRegion_RC = imcrop(img_RC,
            
            % Calculate percentage area overlap with green channel
            overlap_with_GC_mask = and(logical(mask_GC),logical(maskCurrentSynaptosome));
            overlap_with_GC = (sum(overlap_with_GC_mask)/sum(logical(maskCurrentSynaptosome)))*100;
            overlap_perc_GC(j) = overlap_with_GC;
            disp([num2str(overlap_with_GC) '% overlap with green channel.']);

            % Calculate percentage area overlap with blue channel
            overlap_with_BC_mask = and(logical(mask_BC),logical(maskCurrentSynaptosome));
            overlap_with_BC = (sum(overlap_with_BC_mask)/sum(logical(maskCurrentSynaptosome)))*100;
            overlap_perc_BC(j) = overlap_with_BC;
            disp([num2str(overlap_with_BC) '% overlap with blue channel.']);

            % Apply mask to all channels - not using this anymore
            %masked_red_SR_img   = logical(maskCurrentSynaptosome).*double(img_RC);
            masked_red_SR_img   = logical(maskBoundingRegion).*double(img_RC);
            masked_green_SR_img = logical(maskBoundingRegion).*double(img_GC);
            masked_blue_SR_img  = logical(maskBoundingRegion).*double(img_BC);

            % now crop a region with the bounding box for all the three
            % channels and save the images. 
            red_crop{j}   = imcrop(img_RC,boundingregionprops.BoundingBox);
            green_crop{j} = imcrop(img_GC,boundingregionprops.BoundingBox);
            blue_crop{j}  = imcrop(img_BC,boundingregionprops.BoundingBox);
            
            imwrite(flip(flip(red_crop{j},1),2),  fullfile(path_output,'regions',['clusterRC_' num2str(j) '.png']));
            imwrite(flip(flip(green_crop{j},1),2),fullfile(path_output,'regions',['clusterGC_' num2str(j) '.png']));
            imwrite(flip(flip(blue_crop{j},1),2), fullfile(path_output,'regions',['clusterBC_' num2str(j) '.png']));

            % Get weighted overlap for red and green channel
            X = getColocCoefficient(masked_red_SR_img,masked_green_SR_img);
            weighted_overlap_GC(j) = X.mandersCoeff;
            overlap_RonG(j) = X.M1;
            overlap_GonR(j) = X.M2;
            
            % Get weighted overlap for red and blue channel
            X = getColocCoefficient(masked_red_SR_img,masked_blue_SR_img);
            weighted_overlap_BC(j) = X.mandersCoeff;
            overlap_RonB(j) = X.M1;
            overlap_BonR(j) = X.M2; 
            
            % Get weighted overlap for green and blue channel
            X = getColocCoefficient(masked_green_SR_img,masked_blue_SR_img);
            weighted_overlap_GB(j) = X.mandersCoeff;
            overlap_GonB(j) = X.M1;
            overlap_BonG(j) = X.M2; 
            disp(' ');
           
        end

        % Summarize meaurements in a table
        synaptosome_ID = 1:numberOfClusters;
        results = [synaptosome_ID' x_centroid' y_centroid' areas' overlap_perc_GC' overlap_perc_BC' weighted_overlap_GC' weighted_overlap_BC' weighted_overlap_GB' overlap_RonG' overlap_GonR' overlap_RonB' overlap_BonR'];
        results = array2table(results,'VariableNames', {'ID','xCentroid','yCentroid','Area','OverlapWithGreen','OverlapWithBlue','WeightedOverlapWithGreen','WeightedOverlapWithBlue','WeightedOverlapGreenBlue','MandersRG','MandersGR','MandersRB','MandersBR'});
        filename = char(strcat(area_token,'_results.csv'));
        writetable(results,fullfile(path_output, filename));
        fclose(summary_file);
        
        %% Plot results

        % Overlap red and green channel (masks)
        mask_RC_GC = zeros(size(mask_RC,1),size(mask_RC,2),3);
        mask_RC_GC(:,:,1) = mask_RC;
        mask_RC_GC(:,:,2) = mask_GC;
        if show
            figure;
            imshow(flip(flip(mask_RC_GC,2),1));
            title('Red and green channel','FontSize',15)
        end
        if ~essence
            filename = char(strcat(area_token,'_overlap_green.tif'));
            imwrite(flip(flip(mask_RC_GC,2),1), fullfile(path_output, filename));
        end

        % Overlap red and blue channel (masks)
        mask_RC_BC = zeros(size(mask_RC,1),size(mask_RC,2),3);
        mask_RC_BC(:,:,1) = mask_RC;
        mask_RC_BC(:,:,3) = mask_BC;
        if show
            figure;
            imshow(flip(flip(mask_RC_BC,2),1));
            title('Red and blue channel','FontSize',15);
        end
        if ~essence
            name_mask = char(strcat(area_token,'_overlap_blue.tif'));
            imwrite(flip(flip(mask_RC_BC,2),1), fullfile(path_output, name_mask));
        end
    else
        continue
    end
end

toc