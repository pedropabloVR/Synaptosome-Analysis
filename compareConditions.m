%% Synapto-Analysis: compareConditions

% This script reads in the results from the synaptosomeAnalysis.m script
% and further processes the data to obtain a set of results from 
%
% The script outputs:
% 	- Super resolution 3-colour reconstructions of the analyzed regions,
% 	with circles around the detected synaptosomes.
%
% 	- .csv and .mat file with the concatenated analysis results before and
% 	after thresholding based on color channel overlaps.
%
% 	- An estimate of cluster size based on ripley's L(r ) -r function, and
% 	the cluster-to-cluster distances
%
% 	- Figures with plots of the unweighted and weighted overlap
% 	coefficients between the red/green and red/blue channels, the cluster
% 	size estimate from Ripley's K function, the inter-cluster distances.
%
%   - 8-bit reconstructions and localisation files of the cropped regions
%   of interest used for the analysis of Ripley's K function (and all other
%   analysis).

% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated: 6 Feb 2018

clear all
close all
clc
tic

%% File I/O parameters

% File input parameters
repeats = {'A','B'}; %,'Oct-Nov'}; % Name of the two different repeats for each temperature condition

% Path to directories containing results from Synaptosome analysis for different repeats (can be a cell with only one repeat!)
Dir_PHYS      = {fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/Results_4Ca phys'),...
                 fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/Results_4Cb phys')};
Dir_EGTA      = {fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/Results_4Ca egta'),...
                 fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/Results_4Cb egta')};
Dir_EGTAK     = {fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/Results_4Ca egtak'),...
                 fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/Results_4Cb egta')};
        
% Path to raw localisation files for each repeat (no need to change these,
% raw data stays the same) - these are used to calculate the rejection rate. 
dir_phys_raw  = {fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/4Ca_phys'),...
                 fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/4Cb_phys')};
dir_egta_raw  = {fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/4Ca_egta'),...
                 fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/4Cb_egta')};
dir_egtak_raw = {fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/4Ca_egtak'),...
                 fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/4Cb_egtak')};
            
% Output directory
output_dir      = fullfile('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/test data for synapto-analysis/');
format          = 'thunderstorm';

% Reconstructed image parameters
magnification   = 10; % determines final pixel size in super-resolution reconstructions
imgtype         = '_reconstruction_RGB.tif'; % image to mark detected synaptosomes on with circles
radius          = 100; % radius of the circles drawn around detected synaptosomes
colour          = [100 100 100]; % colour of the circles (triplet)

%% Filtering parameters

% Channel overlap filter to detect synaptosomes
filterOverlapGreen = 0;  % 1 to filter on overlap between red and green to detect synaptosomes
filterOverlapBlue  = 0;  % 1 to filter on overlap between red and blue  to detect synaptosomes
minOverlapGreen    = 20;  % minimum overlap between red and green to detect synaptosomes
minOverlapBlue     = 20; % minimum overlap between red and green to detect synaptosomes

% Proximity filter: minimum distance between two synaptosomes
min_dist_between_synaptosomes   = 300; % in nm

% Size filter: very "large" synaptosomes, signaled by large areas of the
% mCLING, are removed
max_Area    = 2000; % max area of mCling

% Ripley's L(r) - r function parameters
r_step      = 10;
R_max       = 1000;
windowsize  = 80; % pixels
pixelsize   = 117; % nm

% Flags for showing results
plot_results_individual_repeats = 1; % 1 to also plot results of individual repeats seperately
show                            = 0; % to show intermediate results
flagprint                       = 1; % set to 1 to save the fig visualization of the results


%% Loop over repeats

clear var results_repeats_pooled

for n = 1:length(repeats)
    
    % Get directories of current repeat
    dir_PHYS    = Dir_PHYS{n};
    dir_EGTA    = Dir_EGTA{n};
    dir_EGTAK   = Dir_EGTAK{n};

    dir_PHYS_raw    = dir_phys_raw{n};
    dir_EGTA_raw    = dir_egta_raw{n};
    dir_EGTAK_raw   = dir_egtak_raw{n};
    
    % Create new output folder
    path_output = fullfile(output_dir,'Results_combined',['Repeat_' repeats{n}]);
    mkdir(path_output);

    save(fullfile(path_output,'parameters.mat'));
    diary(fullfile(path_output,'command_window_output.txt'));

    
    %% Read in results of different ROIs and merge them

    % PHYS ----------------------------------------------------------------

    % Read in all PHYS results and merge in one table
    if ispc
        filelist_PHYS = dir([dir_PHYS '\*results.csv']);
    elseif ismac
        filelist_PHYS = dir([dir_PHYS '/*results.csv']);
    end

    % Initialize
    results_PHYS = [];
    sampleIDs = {};

    for i = 1:length(filelist_PHYS)

        % Get filename and path
        filename = filelist_PHYS(i).name;
        if filename(1) == '.'; continue; end
        filepath = fullfile(dir_PHYS,filename);

        % Get sample ID
        sampleID = strsplit(filename,'_');
        sampleID = sampleID{1};

        % Read in results file
        currentFile = dlmread(filepath,',',1,0);

        % Get number of synaptosomes in current file
        num_synaptosomes = size(currentFile,1);

        % Add current sampleID the same number of times to sampleIDs cell
        clear var sampleID_i
        [sampleID_i{1:num_synaptosomes}] = deal([sampleID '_phys']);
        sampleIDs = [sampleIDs sampleID_i];

        % Concatenate results
        results_PHYS = [results_PHYS; currentFile];
    end

    % Convert concatenated results to table
    results_PHYS = array2table(results_PHYS,'VariableNames',...
        {'synaptosomeID','xCentroid','yCentroid','Area',...
        'OverlapWithGreen','OverlapWithBlue',...
        'WeightedOverlapWithGreen','WeightedOverlapWithBlue'});

    % Add overlapping green and blue area
    results_PHYS.AreaGreen = (results_PHYS.OverlapWithGreen.*results_PHYS.Area)/100;
    results_PHYS.AreaBlue  = (results_PHYS.OverlapWithBlue.*results_PHYS.Area)/100;

    % Add sampleID column in front
    results_PHYS.sampleID  = sampleIDs';
    results_PHYS = results_PHYS(:,[end 1:end-1]);
    
    % EGTA ----------------------------------------------------------------

    % Read in all EGTA results and merge in one table
    if ispc
        filelist_EGTA = dir([dir_EGTA '\*results.csv']);
    elseif ismac
        filelist_EGTA = dir([dir_EGTA '/*results.csv']);
    end

    % Initialize
    results_EGTA = [];
    sampleIDs = {};

    for i = 1:length(filelist_EGTA)

        % Get filename and path
        filename = filelist_EGTA(i).name;
        if filename(1) == '.'; continue; end
        filepath = fullfile(dir_EGTA,filename);

        % Get sample ID
        sampleID = strsplit(filename,'_');
        sampleID = sampleID{1};

        % Read in results file
        currentFile = dlmread(filepath,',',1,0);

        % Get number of synaptosomes in current file
        num_synaptosomes = size(currentFile,1);

        % Add current sampleID the same number of times to sampleIDs cell
        clear var sampleID_i
        [sampleID_i{1:num_synaptosomes}] = deal([sampleID '_egta']);
        sampleIDs = [sampleIDs sampleID_i];

        % Concatenate results
        results_EGTA = [results_EGTA; currentFile];
    end

    % Convert concatenated results to table
    results_EGTA = array2table(results_EGTA,'VariableNames',...
        {'synaptosomeID','xCentroid','yCentroid','Area',...
        'OverlapWithGreen','OverlapWithBlue',...
        'WeightedOverlapWithGreen','WeightedOverlapWithBlue'});

    % Add overlapping green and blue area
    results_EGTA.AreaGreen = (results_EGTA.OverlapWithGreen.*results_EGTA.Area)/100;
    results_EGTA.AreaBlue  = (results_EGTA.OverlapWithBlue.*results_EGTA.Area)/100;

    % Add sampleID column in front
    results_EGTA.sampleID  = sampleIDs';
    results_EGTA = results_EGTA(:,[end 1:end-1]);


    % EGTAK ---------------------------------------------------------------

    % Read in all EGTAK results and merge in one table
    if ispc
        filelist_EGTAK = dir([dir_EGTAK '\*results.csv']);
    elseif ismac
        filelist_EGTAK = dir([dir_EGTAK '/*results.csv']);
    end

    % Initialize
    results_EGTAK = [];
    sampleIDs = {};

    for i = 1:length(filelist_EGTAK)

        % Get filename and path
        filename = filelist_EGTAK(i).name;
        if filename(1) == '.'; continue; end
        filepath = fullfile(dir_EGTAK,filename);

        % Get sample ID
        sampleID = strsplit(filename,'_');
        sampleID = sampleID{1};

        % Read in results file
        currentFile = dlmread(filepath,',',1,0);

        % Get number of synaptosomes in current file
        num_synaptosomes = size(currentFile,1);

        % Add current sampleID the same number of times to sampleIDs cell
        clear var sampleID_i
        [sampleID_i{1:num_synaptosomes}] = deal([sampleID '_egtak']);
        sampleIDs = [sampleIDs sampleID_i];

        % Concatenate results
        results_EGTAK = [results_EGTAK; currentFile];
    end

    % Convert concatenated results to table
    results_EGTAK = array2table(results_EGTAK,'VariableNames',...
        {'synaptosomeID','xCentroid','yCentroid','Area',...
        'OverlapWithGreen','OverlapWithBlue',...
        'WeightedOverlapWithGreen','WeightedOverlapWithBlue'});

    % Add overlapping green and blue area
    results_EGTAK.AreaGreen = (results_EGTAK.OverlapWithGreen.*results_EGTAK.Area)/100;
    results_EGTAK.AreaBlue  = (results_EGTAK.OverlapWithBlue.*results_EGTAK.Area)/100;

    % Add sampleID column in front
    results_EGTAK.sampleID  = sampleIDs';
    results_EGTAK = results_EGTAK(:,[end 1:end-1]);


    %% Merge tables of the three conditions before filtering

    results_combined = [results_PHYS; results_EGTA; results_EGTAK];
    
    clear var phys egta egtak
    [phys{1:size(results_PHYS,1)}] = deal('phys');
    [egta{1:size(results_EGTA,1)}] = deal('egta');
    [egtak{1:size(results_EGTAK,1)}] = deal('egtak');
    condition_column = [phys egta egtak];

    results_combined.condition = condition_column';
    results_combined = results_combined(:,[end 1:end-1]);

    path_results_combined = fullfile(path_output,'results_combined_before_overlap_threshold.mat');
    save(path_results_combined,'results_combined');

    path_results_combined = fullfile(path_output,'results_combined_before_overlap_threshold.txt');
    writetable(results_combined,path_results_combined,'Delimiter','\t');


    %% Overlap filter: Detect synaptosomes by filtering on r/g and r/b overlap
    
    % Many of the blobs in the red channel for which overlaps are calculated
    % are not actually synaptosomes. If they don'toverlap with blobs in the
    % blue channel, they probably aren't synaptosomes. So, we filter out all
    % results for red blobs that don't show significant overlap with the green
    % channel.

    synapto_number_prefilter_PHYS = size(results_PHYS);
    synapto_number_prefilter_EGTA = size(results_EGTA);
    synapto_number_prefilter_EGTAK = size(results_EGTAK);

    if filterOverlapGreen
        results_PHYS  = results_PHYS(results_PHYS.OverlapWithGreen   > minOverlapGreen,:);
        results_EGTA  = results_EGTA(results_EGTA.OverlapWithGreen   > minOverlapGreen,:);
        results_EGTAK = results_EGTAK(results_EGTAK.OverlapWithGreen > minOverlapGreen,:);
    end 
    % The elseif in the previous loop was never executed as long as
    % filterOverlapGreen was true
    if filterOverlapBlue

        results_PHYS  = results_PHYS(results_PHYS.OverlapWithBlue   > minOverlapBlue,:);
        results_EGTA  = results_EGTA(results_EGTA.OverlapWithBlue   > minOverlapBlue,:);
        results_EGTAK = results_EGTAK(results_EGTAK.OverlapWithBlue > minOverlapBlue,:);
    end

    synapto_number_post_overlapFilter_PHYS = size(results_PHYS);
    synapto_number_post_overlapFilter_EGTA = size(results_EGTA);
    synapto_number_post_overlapFilter_EGTAK = size(results_EGTAK);

    disp(['Synaptosomes PHYS pre filter: ' num2str(synapto_number_prefilter_PHYS(1))]);
    disp(['Synaptosomes PHYS post filter: ' num2str(synapto_number_post_overlapFilter_PHYS(1))]);
    disp(['Synaptosomes EGTA pre filter: ' num2str(synapto_number_prefilter_EGTA(1))]);
    disp(['Synaptosomes EGTA post filter: ' num2str(synapto_number_post_overlapFilter_EGTA(1))]);
    disp(['Synaptosomes EGTAK pre filter: ' num2str(synapto_number_prefilter_EGTAK(1))]);
    disp(['Synaptosomes EGTAK post filter: ' num2str(synapto_number_post_overlapFilter_EGTAK(1))]);


    %% Proximity filter: Remove synaptosomes that are too close together

    % Make a copy of the results before the proximity filtering
    % (so only the overlap-filtered results)
    results_PHYS_before_proximity_filter = results_PHYS;
    results_EGTA_before_proximity_filter = results_EGTA;
    results_EGTAK_before_proximity_filter = results_EGTAK;

    clear var filtered_results_PHYS
    clear var filtered_results_EGTA
    clear var filtered_results_EGTAK
    
    % PHYS condition ------------------------------------------------------
    
    % Loop over sample_ID in the results
    sample_IDs = unique(results_PHYS.sampleID);
    for i=1:size(sample_IDs,1)

        % Get synaptosomes from current sample
        sample_ID_i = sample_IDs{i};
        results_sampleID = results_PHYS(strcmp(results_PHYS.sampleID,sample_ID_i),:);
        synaptosome_number_pre_filter = size(results_sampleID.synaptosomeID);

        % Get indeces of synaptosomes in this sample that are too close together
        points = [results_sampleID.xCentroid results_sampleID.yCentroid];
        pairwise_dist_matrix = pdist2(points,points);
        to_remove = (pairwise_dist_matrix < min_dist_between_synaptosomes) - eye(size(points,1));
        index_to_remove = sum(to_remove) > 0;
        disp(['Condition and sample ID: ' sample_ID_i]);
        disp(['Synaptosomes before proximity filter: ' num2str(synaptosome_number_pre_filter(1))]);
        disp(['To remove: ' num2str(sum(index_to_remove))]);

        % Remove the synaptosomes that are too close together
        filtered_results_sampleID = results_sampleID(index_to_remove == 0,:);
        if ~exist('filtered_results_PHYS')
            filtered_results_PHYS = filtered_results_sampleID;
        else
            filtered_results_PHYS = [filtered_results_PHYS; filtered_results_sampleID];
        end
    end

    results_PHYS = filtered_results_PHYS;
    

    % EGTA condition ------------------------------------------------------

    % Loop over sample_ID in the results
    sample_IDs = unique(results_EGTA.sampleID);
    for i=1:size(sample_IDs,1)

        % Get synaptosomes from current sample
        sample_ID_i = sample_IDs{i};
        results_sampleID = results_EGTA(strcmp(results_EGTA.sampleID,sample_ID_i),:);
        synaptosome_number_pre_filter = size(results_sampleID.synaptosomeID);

        % Get indeces of synaptosomes in this sample that are too close together
        points = [results_sampleID.xCentroid results_sampleID.yCentroid];
        pairwise_dist_matrix = pdist2(points,points);
        to_remove = (pairwise_dist_matrix < min_dist_between_synaptosomes) - eye(size(points,1));
        index_to_remove = sum(to_remove) > 0;
        disp(['Condition and sample ID: ' sample_ID_i]);
        disp(['Synaptosomes before proximity filter: ' num2str(synaptosome_number_pre_filter(1))]);
        disp(['To remove: ' num2str(sum(index_to_remove))]);

        % Remove the synaptosomes that are too close together
        filtered_results_sampleID = results_sampleID(index_to_remove == 0,:);
        if ~exist('filtered_results_EGTA')
            filtered_results_EGTA = filtered_results_sampleID;
        else
            filtered_results_EGTA = [filtered_results_EGTA; filtered_results_sampleID];
        end
    end

    results_EGTA = filtered_results_EGTA;

    
    % EGTAK condition ------------------------------------------------------

    % Loop over sample_ID in the results
    sample_IDs = unique(results_EGTAK.sampleID);
    for i=1:size(sample_IDs,1)

        % Get synaptosomes from current sample
        sample_ID_i = sample_IDs{i};
        results_sampleID = results_EGTAK(strcmp(results_EGTAK.sampleID,sample_ID_i),:);
        synaptosome_number_pre_filter = size(results_sampleID.synaptosomeID);

        % Get indeces of synaptosomes in this sample that are too close together
        points = [results_sampleID.xCentroid results_sampleID.yCentroid];
        pairwise_dist_matrix = pdist2(points,points);
        to_remove = (pairwise_dist_matrix < min_dist_between_synaptosomes) - eye(size(points,1));
        index_to_remove = sum(to_remove) > 0;
        disp(['Condition and sample ID: ' sample_ID_i]);
        disp(['Synaptosomes before proximity filter: ' num2str(synaptosome_number_pre_filter(1))]);
        disp(['To remove: ' num2str(sum(index_to_remove))]);

        % Remove the synaptosomes that are too close together
        filtered_results_sampleID = results_sampleID(index_to_remove == 0,:);
        if ~exist('filtered_results_EGTAK')
            filtered_results_EGTAK = filtered_results_sampleID;
        else
            filtered_results_EGTAK = [filtered_results_EGTAK; filtered_results_sampleID];
        end
    end

    results_EGTAK = filtered_results_EGTAK;

    % -------------------------------------------------------------------------

    synapto_number_post_proximityFilter_PHYS = size(results_PHYS);
    synapto_number_post_proximityFilter_EGTA = size(results_EGTA);
    synapto_number_post_proximityFilter_EGTAK = size(results_EGTAK);

    disp(['Synaptosomes PHYS post proximity filter: '  num2str(synapto_number_post_proximityFilter_PHYS(1))]);
    disp(['Synaptosomes EGTA post proximity filter: '  num2str(synapto_number_post_proximityFilter_EGTA(1))]);
    disp(['Synaptosomes EGTAK post proximity filter: ' num2str(synapto_number_post_proximityFilter_EGTAK(1))]);


    %% Filter on area of mCling
    
    results_PHYS  = results_PHYS(results_PHYS.Area < max_Area,:);
    results_EGTA  = results_EGTA(results_EGTA.Area < max_Area,:);
    results_EGTAK = results_EGTAK(results_EGTAK.Area < max_Area,:);
    
    %% Add column to record number of localisations per region analyzed in Ripley's K
    
    num_regions_PHYS = size(results_PHYS,1); 
    num_regions_EGTA = size(results_EGTA,1); 
    num_regions_EGTAK = size(results_EGTAK,1); 
    
    results_PHYS.numLocs_RC = zeros(num_regions_PHYS,1);
    results_EGTA.numLocs_RC = zeros(num_regions_EGTA,1);
    results_EGTAK.numLocs_RC = zeros(num_regions_EGTAK,1);

    results_PHYS.numLocs_GC = zeros(num_regions_PHYS,1);
    results_EGTA.numLocs_GC = zeros(num_regions_EGTA,1);
    results_EGTAK.numLocs_GC = zeros(num_regions_EGTAK,1);
    
    results_PHYS.numLocs_BC = zeros(num_regions_PHYS,1);
    results_EGTA.numLocs_BC = zeros(num_regions_EGTA,1);
    results_EGTAK.numLocs_BC = zeros(num_regions_EGTAK,1);
    
    
    %% Ripley's K analysis

    % Create new subfolder in output folder
    path_output_ripley = fullfile(path_output,'ripley');
    mkdir(path_output_ripley);
    mkdir(fullfile(path_output_ripley,'phys'));
    mkdir(fullfile(path_output_ripley,'egta'));
    mkdir(fullfile(path_output_ripley,'egtak'));

    % PHYS condition ------------------------------------------------------

    condition = 'phys';

    % Initialize
    ripley_RR = [];
    ripley_GG = [];
    ripley_BB = [];
    ripley_RG = [];
    ripley_RB = [];
    ripley_GB = [];
    
    numLocs_RC_PHYS = [];
    numLocs_GC_PHYS = [];
    numLocs_BC_PHYS = [];

    rejectionRate_RC_PHYS = [];
    rejectionRate_GC_PHYS = [];
    rejectionRate_BC_PHYS = [];
    
    % Create subfolders
    mkdir(fullfile(path_output_ripley,condition,'RC'));
    mkdir(fullfile(path_output_ripley,condition,'GC'));
    mkdir(fullfile(path_output_ripley,condition,'BC'));
    mkdir(fullfile(path_output_ripley,condition,'images'));

    % Get sample IDs and loop over them
    sampleIDs = unique(results_PHYS.sampleID);
    for i = 1:length(sampleIDs)

        % Load localisation files of sample filtered by fitting parameters.
        area_token = strsplit(char(sampleIDs(i)),'_');
        area_token = area_token{1};
        locs_RC = readLocFile(fullfile(dir_PHYS,strcat(area_token,'_locs_RC_filtered.csv')),format);
        locs_GC = readLocFile(fullfile(dir_PHYS,strcat(area_token,'_locs_GC_filtered.csv')),format);
        locs_BC = readLocFile(fullfile(dir_PHYS,strcat(area_token,'_locs_BC_filtered.csv')),format);
        
        % Load localisation files prior to filtering
        
        locs_RC_prefilter = readLocFile(fullfile(dir_PHYS_raw,strcat(area_token,'_647.csv')),format);
        locs_GC_prefilter = readLocFile(fullfile(dir_PHYS_raw,strcat(area_token,'_561_reg.csv')),format);
        locs_BC_prefilter = readLocFile(fullfile(dir_PHYS_raw,strcat(area_token,'_488_reg.csv')),format);

        % Estimate Fov
        Fov = estimateFov([locs_RC.x; locs_GC.x; locs_GC.x; locs_RC.y; locs_BC.y; locs_BC.y], pixelsize);

        % Get number of synaptosomes and coordinates of their centroids
        results_PHYS_i = results_PHYS(strcmp(results_PHYS.sampleID,sampleIDs(i)),:);
        synaptosomeID = results_PHYS_i.synaptosomeID;
        num_synaptosomes = size(results_PHYS_i,1);
        X_c = results_PHYS_i.xCentroid;
        Y_c = results_PHYS_i.yCentroid;

        % Initialize
        r_hist = 0:r_step:R_max;
        ripley_RR_i = zeros(length(r_hist),num_synaptosomes);
        ripley_GG_i = zeros(length(r_hist),num_synaptosomes);
        ripley_BB_i = zeros(length(r_hist),num_synaptosomes);
        ripley_RG_i = zeros(length(r_hist),num_synaptosomes);
        ripley_RB_i = zeros(length(r_hist),num_synaptosomes);
        ripley_GB_i = zeros(length(r_hist),num_synaptosomes);    
        
        num_Locs_RC_i = zeros(num_synaptosomes,1);
        num_Locs_GC_i = zeros(num_synaptosomes,1);
        num_Locs_BC_i = zeros(num_synaptosomes,1);
        
        rejectionRate_RC_PHYS_i = zeros(num_synaptosomes,1);
        rejectionRate_GC_PHYS_i = zeros(num_synaptosomes,1);
        rejectionRate_BC_PHYS_i = zeros(num_synaptosomes,1);
        
        % Loop over all synaptosomes
        for j = 1:num_synaptosomes

            % Get coordinates of synaptosome
            x_c = X_c(j)*magnification;
            y_c = Y_c(j)*magnification;

            % Crop locfiles to region around synaptosome
            locs_RC_cropped = cropLocsAroundCentroid(locs_RC,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_GC_cropped = cropLocsAroundCentroid(locs_GC,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_BC_cropped = cropLocsAroundCentroid(locs_BC,x_c,y_c,(windowsize/magnification)*pixelsize);
            
            locs_RC_cropped_prefilter = cropLocsAroundCentroid(locs_RC_prefilter,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_GC_cropped_prefilter = cropLocsAroundCentroid(locs_GC_prefilter,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_BC_cropped_prefilter = cropLocsAroundCentroid(locs_BC_prefilter,x_c,y_c,(windowsize/magnification)*pixelsize);
           
            % rejection rate 
            rejectionRate_RC_PHYS_i(j,1) = (size(locs_RC_cropped_prefilter,1) - size(locs_RC_cropped,1))/size(locs_RC_cropped_prefilter,1) * 100;
            rejectionRate_GC_PHYS_i(j,1) = (size(locs_GC_cropped_prefilter,1) - size(locs_GC_cropped,1))/size(locs_GC_cropped_prefilter,1) * 100;
            rejectionRate_BC_PHYS_i(j,1) = (size(locs_BC_cropped_prefilter,1) - size(locs_BC_cropped,1))/size(locs_BC_cropped_prefilter,1) * 100;
            
            % Write away cropped localisation files
            writeLocFile(locs_RC_cropped,fullfile(path_output_ripley,condition,'RC',strcat(area_token,'_',condition,'_',num2str(synaptosomeID(j)),'.csv')),format)
            writeLocFile(locs_GC_cropped,fullfile(path_output_ripley,condition,'GC',strcat(area_token,'_',condition,'_',num2str(synaptosomeID(j)),'.csv')),format)
            writeLocFile(locs_BC_cropped,fullfile(path_output_ripley,condition,'BC',strcat(area_token,'_',condition,'_',num2str(synaptosomeID(j)),'.csv')),format)

            X_RC = locs_RC_cropped.x;
            Y_RC = locs_RC_cropped.y;
            X_GC = locs_GC_cropped.x;
            Y_GC = locs_GC_cropped.y;
            X_BC = locs_BC_cropped.x;
            Y_BC = locs_BC_cropped.y;
            
            num_Locs_RC_i(j,1) = size(X_RC,1);
            num_Locs_GC_i(j,1) = size(X_GC,1);
            num_Locs_BC_i(j,1) = size(X_BC,1);

            % Get images of cropped regions and write them away
            % These images will be used later in the synapto_display.m script
            sigma = 150/magnification;
            intensities_RC = locs_RC_cropped.intensity;
            intensities_GC = locs_GC_cropped.intensity;
            intensities_BC = locs_BC_cropped.intensity;

            img_RC = generateImage(X_RC-min(X_RC)+1, Y_RC-min(Y_RC)+1, sigma, intensities_RC, windowsize/magnification, pixelsize, magnification);
            img_GC = generateImage(X_GC-min(X_GC)+1, Y_GC-min(Y_GC)+1, sigma, intensities_GC, windowsize/magnification, pixelsize, magnification);
            img_BC = generateImage(X_BC-min(X_BC)+1, Y_BC-min(Y_BC)+1, sigma, intensities_BC, windowsize/magnification, pixelsize, magnification);

            img_name_RC = strcat(char(sampleIDs(i)),'_RC_synaptosome_',num2str(synaptosomeID(j)),'.png');
            img_name_GC = strcat(char(sampleIDs(i)),'_GC_synaptosome_',num2str(synaptosomeID(j)),'.png');
            img_name_BC = strcat(char(sampleIDs(i)),'_BC_synaptosome_',num2str(synaptosomeID(j)),'.png');

            imwrite(flip(flip(img_RC,1),2),fullfile(path_output_ripley,condition,'images',img_name_RC));
            imwrite(flip(flip(img_GC,1),2),fullfile(path_output_ripley,condition,'images',img_name_GC));
            imwrite(flip(flip(img_BC,1),2),fullfile(path_output_ripley,condition,'images',img_name_BC));

            % Ripley's K single colour analysis
            [~,~,K_R,~,~] = CalcRipleyK(X_RC,Y_RC, X_RC,Y_RC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_G,~,~] = CalcRipleyK(X_GC,Y_GC, X_GC,Y_GC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_B,~,~] = CalcRipleyK(X_BC,Y_BC, X_BC,Y_BC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            ripley_RR_i(:,j) = sqrt(K_R'/pi) - r_hist;
            ripley_GG_i(:,j) = sqrt(K_G'/pi) - r_hist;
            ripley_BB_i(:,j) = sqrt(K_B'/pi) - r_hist;

            % Ripley's K two colour analysis
            [~,~,K_RG,~,~] = CalcRipleyK(X_RC,Y_RC, X_GC,Y_GC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_RB,~,~] = CalcRipleyK(X_RC,Y_RC, X_BC,Y_BC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_GB,~,~] = CalcRipleyK(X_GC,Y_GC, X_BC,Y_BC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            ripley_RG_i(:,j) = sqrt(K_RG'/pi) - r_hist;
            ripley_RB_i(:,j) = sqrt(K_RB'/pi) - r_hist;
            ripley_GB_i(:,j) = sqrt(K_GB'/pi) - r_hist;
        end
        % Concatenate results ripley with previous samples
        ripley_RR = [ripley_RR ripley_RR_i];
        ripley_GG = [ripley_GG ripley_GG_i];
        ripley_BB = [ripley_BB ripley_BB_i];
        ripley_RG = [ripley_RG ripley_RG_i];
        ripley_RB = [ripley_RB ripley_RB_i];
        ripley_GB = [ripley_GB ripley_GB_i];
        
        numLocs_RC_PHYS = [numLocs_RC_PHYS num_Locs_RC_i'];
        numLocs_GC_PHYS = [numLocs_GC_PHYS num_Locs_GC_i'];
        numLocs_BC_PHYS = [numLocs_BC_PHYS num_Locs_BC_i'];
        
        rejectionRate_RC_PHYS = [rejectionRate_RC_PHYS rejectionRate_RC_PHYS_i'];
        rejectionRate_GC_PHYS = [rejectionRate_GC_PHYS rejectionRate_GC_PHYS_i'];
        rejectionRate_BC_PHYS = [rejectionRate_BC_PHYS rejectionRate_BC_PHYS_i'];
        
        
    end

    H_all_phys_RC = ripley_RR;
    H_all_phys_GC = ripley_GG;
    H_all_phys_BC = ripley_BB;
    H_all_phys_RG = ripley_RG;
    H_all_phys_RB = ripley_RB;
    H_all_phys_GB = ripley_GB;

    results_PHYS.numLocs_RC = numLocs_RC_PHYS';
    results_PHYS.numLocs_GC = numLocs_GC_PHYS';
    results_PHYS.numLocs_BC = numLocs_BC_PHYS';
    
    results_PHYS.rejectionRate_RC = rejectionRate_RC_PHYS';
    results_PHYS.rejectionRate_GC = rejectionRate_GC_PHYS';
    results_PHYS.rejectionRate_BC = rejectionRate_BC_PHYS';
        
    
    % Write away Ripley's functions before removing NaN columns (to retain
    % information about from which synaptosome the curve came)
    path_ripley_phys = fullfile(path_output_ripley,'H_phys.mat');
    save(path_ripley_phys,'ripley_RR','ripley_GG','ripley_BB',...
                          'ripley_RG','ripley_RB','ripley_GB');

    % clear unfiltered localisation tables to free memory
    clear locs_RC_prefilter locs_GC_prefilter locs_BC_prefilter;
    

    % EGTA condition ------------------------------------------------------

    condition = 'egta';

    % Initialize
    ripley_RR = [];
    ripley_GG = [];
    ripley_BB = [];
    ripley_RG = [];
    ripley_RB = [];
    ripley_GB = [];
    
    numLocs_RC_EGTA = [];
    numLocs_GC_EGTA = [];
    numLocs_BC_EGTA = [];

    rejectionRate_RC_EGTA = [];
    rejectionRate_GC_EGTA = [];
    rejectionRate_BC_EGTA = [];
    
    % Create subfolders
    mkdir(fullfile(path_output_ripley,condition,'RC'));
    mkdir(fullfile(path_output_ripley,condition,'GC'));
    mkdir(fullfile(path_output_ripley,condition,'BC'));
    mkdir(fullfile(path_output_ripley,condition,'images'));

    % Get sample IDs and loop over them
    sampleIDs = unique(results_EGTA.sampleID);
    for i = 1:length(sampleIDs)

        % Load localisation files of sample
        area_token = strsplit(char(sampleIDs(i)),'_');
        area_token = area_token{1};
        locs_RC = readLocFile(fullfile(dir_EGTA,strcat(area_token,'_locs_RC_filtered.csv')),format);
        locs_GC = readLocFile(fullfile(dir_EGTA,strcat(area_token,'_locs_GC_filtered.csv')),format);
        locs_BC = readLocFile(fullfile(dir_EGTA,strcat(area_token,'_locs_BC_filtered.csv')),format);

        locs_RC_prefilter = readLocFile(fullfile(dir_EGTA_raw,strcat(area_token,'_647.csv')),format);
        locs_GC_prefilter = readLocFile(fullfile(dir_EGTA_raw,strcat(area_token,'_561_reg.csv')),format);
        locs_BC_prefilter = readLocFile(fullfile(dir_EGTA_raw,strcat(area_token,'_488_reg.csv')),format);

        % Estimate Fov
        Fov = estimateFov([locs_RC.x; locs_GC.x; locs_GC.x; locs_RC.y; locs_BC.y; locs_BC.y], pixelsize);

        % Get number of synaptosomes and coordinates of their centroids
        results_EGTA_i = results_EGTA(strcmp(results_EGTA.sampleID,sampleIDs(i)),:);
        synaptosomeID = results_EGTA_i.synaptosomeID;
        num_synaptosomes = size(results_EGTA_i,1);
        X_c = results_EGTA_i.xCentroid;
        Y_c = results_EGTA_i.yCentroid;

        % Initialize
        r_hist = 0:r_step:R_max;
        ripley_RR_i = zeros(length(r_hist),num_synaptosomes);
        ripley_GG_i = zeros(length(r_hist),num_synaptosomes);
        ripley_BB_i = zeros(length(r_hist),num_synaptosomes);
        ripley_RG_i = zeros(length(r_hist),num_synaptosomes);
        ripley_RB_i = zeros(length(r_hist),num_synaptosomes);
        ripley_GB_i = zeros(length(r_hist),num_synaptosomes);    

        num_Locs_RC_i = zeros(num_synaptosomes,1);
        num_Locs_GC_i = zeros(num_synaptosomes,1);
        num_Locs_BC_i = zeros(num_synaptosomes,1);
                
        % Loop over all synaptosomes
        for j = 1:num_synaptosomes

            % Get coordinates of synaptosome
            x_c = X_c(j)*magnification;
            y_c = Y_c(j)*magnification;

            % Crop locfiles to region around synaptosome
            locs_RC_cropped = cropLocsAroundCentroid(locs_RC,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_GC_cropped = cropLocsAroundCentroid(locs_GC,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_BC_cropped = cropLocsAroundCentroid(locs_BC,x_c,y_c,(windowsize/magnification)*pixelsize);
            
            % Crop unfiltered localisation files
            locs_RC_cropped_prefilter = cropLocsAroundCentroid(locs_RC_prefilter,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_GC_cropped_prefilter = cropLocsAroundCentroid(locs_GC_prefilter,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_BC_cropped_prefilter = cropLocsAroundCentroid(locs_BC_prefilter,x_c,y_c,(windowsize/magnification)*pixelsize);
           
            % rejection rate 
            rejectionRate_RC_EGTA_i(j,1) = (size(locs_RC_cropped_prefilter,1) - size(locs_RC_cropped,1))/size(locs_RC_cropped_prefilter,1) * 100;
            rejectionRate_GC_EGTA_i(j,1) = (size(locs_GC_cropped_prefilter,1) - size(locs_GC_cropped,1))/size(locs_GC_cropped_prefilter,1) * 100;
            rejectionRate_BC_EGTA_i(j,1) = (size(locs_BC_cropped_prefilter,1) - size(locs_BC_cropped,1))/size(locs_BC_cropped_prefilter,1) * 100;
            
            % Write away cropped localisation files
            writeLocFile(locs_RC_cropped,fullfile(path_output_ripley,condition,'RC',strcat(area_token,'_',condition,'_',num2str(synaptosomeID(j)),'.csv')),format)
            writeLocFile(locs_GC_cropped,fullfile(path_output_ripley,condition,'GC',strcat(area_token,'_',condition,'_',num2str(synaptosomeID(j)),'.csv')),format)
            writeLocFile(locs_BC_cropped,fullfile(path_output_ripley,condition,'BC',strcat(area_token,'_',condition,'_',num2str(synaptosomeID(j)),'.csv')),format)

            X_RC = locs_RC_cropped.x;
            Y_RC = locs_RC_cropped.y;
            X_GC = locs_GC_cropped.x;
            Y_GC = locs_GC_cropped.y;
            X_BC = locs_BC_cropped.x;
            Y_BC = locs_BC_cropped.y;

            num_Locs_RC_i(j,1) = size(X_RC,1);
            num_Locs_GC_i(j,1) = size(X_GC,1);
            num_Locs_BC_i(j,1) = size(X_BC,1);
            
            % Get images of cropped regions and write them away
            sigma = 150/magnification;
            intensities_RC = locs_RC_cropped.intensity;
            intensities_GC = locs_GC_cropped.intensity;
            intensities_BC = locs_BC_cropped.intensity;

            img_RC = generateImage(X_RC-min(X_RC)+1, Y_RC-min(Y_RC)+1, sigma, intensities_RC, windowsize/magnification, pixelsize, magnification);
            img_GC = generateImage(X_GC-min(X_GC)+1, Y_GC-min(Y_GC)+1, sigma, intensities_GC, windowsize/magnification, pixelsize, magnification);
            img_BC = generateImage(X_BC-min(X_BC)+1, Y_BC-min(Y_BC)+1, sigma, intensities_BC, windowsize/magnification, pixelsize, magnification);

            img_name_RC = strcat(char(sampleIDs(i)),'_RC_synaptosome_',num2str(synaptosomeID(j)),'.png');
            img_name_GC = strcat(char(sampleIDs(i)),'_GC_synaptosome_',num2str(synaptosomeID(j)),'.png');
            img_name_BC = strcat(char(sampleIDs(i)),'_BC_synaptosome_',num2str(synaptosomeID(j)),'.png');

            imwrite(flip(flip(img_RC,1),2),fullfile(path_output_ripley,condition,'images',img_name_RC));
            imwrite(flip(flip(img_GC,1),2),fullfile(path_output_ripley,condition,'images',img_name_GC));
            imwrite(flip(flip(img_BC,1),2),fullfile(path_output_ripley,condition,'images',img_name_BC));

            % Ripley's K single colour analysis
            [~,~,K_R,~,~] = CalcRipleyK(X_RC,Y_RC, X_RC,Y_RC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_G,~,~] = CalcRipleyK(X_GC,Y_GC, X_GC,Y_GC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_B,~,~] = CalcRipleyK(X_BC,Y_BC, X_BC,Y_BC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            ripley_RR_i(:,j) = sqrt(K_R'/pi) - r_hist;
            ripley_GG_i(:,j) = sqrt(K_G'/pi) - r_hist;
            ripley_BB_i(:,j) = sqrt(K_B'/pi) - r_hist;

            % Ripley's K two colour analysis
            [~,~,K_RG,~,~] = CalcRipleyK(X_RC,Y_RC, X_GC,Y_GC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_RB,~,~] = CalcRipleyK(X_RC,Y_RC, X_BC,Y_BC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_GB,~,~] = CalcRipleyK(X_GC,Y_GC, X_BC,Y_BC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            ripley_RG_i(:,j) = sqrt(K_RG'/pi) - r_hist;
            ripley_RB_i(:,j) = sqrt(K_RB'/pi) - r_hist;
            ripley_GB_i(:,j) = sqrt(K_GB'/pi) - r_hist;
        end
        % Concatenate results ripley with previous samples
        ripley_RR = [ripley_RR ripley_RR_i];
        ripley_GG = [ripley_GG ripley_GG_i];
        ripley_BB = [ripley_BB ripley_BB_i];
        ripley_RG = [ripley_RG ripley_RG_i];
        ripley_RB = [ripley_RB ripley_RB_i];
        ripley_GB = [ripley_GB ripley_GB_i];
        
        numLocs_RC_EGTA = [numLocs_RC_EGTA num_Locs_RC_i'];
        numLocs_GC_EGTA = [numLocs_GC_EGTA num_Locs_GC_i'];
        numLocs_BC_EGTA = [numLocs_BC_EGTA num_Locs_BC_i'];
        
        rejectionRate_RC_EGTA = [rejectionRate_RC_EGTA rejectionRate_RC_EGTA_i'];
        rejectionRate_GC_EGTA = [rejectionRate_GC_EGTA rejectionRate_GC_EGTA_i'];
        rejectionRate_BC_EGTA = [rejectionRate_BC_EGTA rejectionRate_BC_EGTA_i'];
        
    end

    H_all_egta_RC = ripley_RR;
    H_all_egta_GC = ripley_GG;
    H_all_egta_BC = ripley_BB;
    H_all_egta_RG = ripley_RG;
    H_all_egta_RB = ripley_RB;
    H_all_egta_GB = ripley_GB;
    
    results_EGTA.numLocs_RC = numLocs_RC_EGTA';
    results_EGTA.numLocs_GC = numLocs_GC_EGTA';
    results_EGTA.numLocs_BC = numLocs_BC_EGTA';

    results_EGTA.rejectionRate_RC = rejectionRate_RC_EGTA';
    results_EGTA.rejectionRate_GC = rejectionRate_GC_EGTA';
    results_EGTA.rejectionRate_BC = rejectionRate_BC_EGTA';
        
    % Write away Ripley's functions before removing NaN columns (to retain
    % information about from which synaptosome the curve came)
    path_ripley_egta = fullfile(path_output_ripley,'H_egta.mat');
    save(path_ripley_egta,'ripley_RR','ripley_GG','ripley_BB',...
                          'ripley_RG','ripley_RB','ripley_GB');
    
    % clear unfiltered localisation tables to free memory
    clear locs_RC_prefilter locs_GC_prefilter locs_BC_prefilter;
    

    % EGTAK condition -----------------------------------------------------

    condition = 'egtak';

    % Initialize
    ripley_RR = [];
    ripley_GG = [];
    ripley_BB = [];
    ripley_RG = [];
    ripley_RB = [];
    ripley_GB = [];
    
    numLocs_RC_EGTAK = [];
    numLocs_GC_EGTAK = [];
    numLocs_BC_EGTAK = [];

    rejectionRate_RC_EGTAK = [];
    rejectionRate_GC_EGTAK = [];
    rejectionRate_BC_EGTAK = [];
    
    % Create subfolders
    mkdir(fullfile(path_output_ripley,condition,'RC'));
    mkdir(fullfile(path_output_ripley,condition,'GC'));
    mkdir(fullfile(path_output_ripley,condition,'BC'));
    mkdir(fullfile(path_output_ripley,condition,'images'));

    % Get sample IDs and loop over them
    sampleIDs = unique(results_EGTAK.sampleID);
    for i = 1:length(sampleIDs)

        % Load localisation files of sample
        area_token = strsplit(char(sampleIDs(i)),'_');
        area_token = area_token{1};
        locs_RC = readLocFile(fullfile(dir_EGTAK,strcat(area_token,'_locs_RC_filtered.csv')),format);
        locs_GC = readLocFile(fullfile(dir_EGTAK,strcat(area_token,'_locs_GC_filtered.csv')),format);
        locs_BC = readLocFile(fullfile(dir_EGTAK,strcat(area_token,'_locs_BC_filtered.csv')),format);

        locs_RC_prefilter = readLocFile(fullfile(dir_EGTAK_raw,strcat(area_token,'_647.csv')),format);
        locs_GC_prefilter = readLocFile(fullfile(dir_EGTAK_raw,strcat(area_token,'_561_reg.csv')),format);
        locs_BC_prefilter = readLocFile(fullfile(dir_EGTAK_raw,strcat(area_token,'_488_reg.csv')),format);

        % Estimate Fov
        Fov = estimateFov([locs_RC.x; locs_GC.x; locs_GC.x; locs_RC.y; locs_BC.y; locs_BC.y], pixelsize);

        % Get number of synaptosomes and coordinates of their centroids
        results_EGTAK_i = results_EGTAK(strcmp(results_EGTAK.sampleID,sampleIDs(i)),:);
        synaptosomeID = results_EGTAK_i.synaptosomeID;
        num_synaptosomes = size(results_EGTAK_i,1);
        X_c = results_EGTAK_i.xCentroid;
        Y_c = results_EGTAK_i.yCentroid;

        % Initialize
        r_hist = 0:r_step:R_max;
        ripley_RR_i = zeros(length(r_hist),num_synaptosomes);
        ripley_GG_i = zeros(length(r_hist),num_synaptosomes);
        ripley_BB_i = zeros(length(r_hist),num_synaptosomes);
        ripley_RG_i = zeros(length(r_hist),num_synaptosomes);
        ripley_RB_i = zeros(length(r_hist),num_synaptosomes);
        ripley_GB_i = zeros(length(r_hist),num_synaptosomes);    

        num_Locs_RC_i = zeros(num_synaptosomes,1);
        num_Locs_GC_i = zeros(num_synaptosomes,1);
        num_Locs_BC_i = zeros(num_synaptosomes,1);
        
        % Loop over all synaptosomes
        for j = 1:num_synaptosomes

            % Get coordinates of synaptosome
            x_c = X_c(j)*magnification;
            y_c = Y_c(j)*magnification;

            % Crop locfiles to region around synaptosome
            locs_RC_cropped = cropLocsAroundCentroid(locs_RC,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_GC_cropped = cropLocsAroundCentroid(locs_GC,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_BC_cropped = cropLocsAroundCentroid(locs_BC,x_c,y_c,(windowsize/magnification)*pixelsize);
            
            % Crop unfiltered localisation files
            locs_RC_cropped_prefilter = cropLocsAroundCentroid(locs_RC_prefilter,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_GC_cropped_prefilter = cropLocsAroundCentroid(locs_GC_prefilter,x_c,y_c,(windowsize/magnification)*pixelsize);
            locs_BC_cropped_prefilter = cropLocsAroundCentroid(locs_BC_prefilter,x_c,y_c,(windowsize/magnification)*pixelsize);
           
            % rejection rate 
            rejectionRate_RC_EGTAK_i(j,1) = (size(locs_RC_cropped_prefilter,1) - size(locs_RC_cropped,1))/size(locs_RC_cropped_prefilter,1) * 100;
            rejectionRate_GC_EGTAK_i(j,1) = (size(locs_GC_cropped_prefilter,1) - size(locs_GC_cropped,1))/size(locs_GC_cropped_prefilter,1) * 100;
            rejectionRate_BC_EGTAK_i(j,1) = (size(locs_BC_cropped_prefilter,1) - size(locs_BC_cropped,1))/size(locs_BC_cropped_prefilter,1) * 100;
            
            % Write away cropped localisation files
            writeLocFile(locs_RC_cropped,fullfile(path_output_ripley,condition,'RC',strcat(area_token,'_',condition,'_',num2str(synaptosomeID(j)),'.csv')),format)
            writeLocFile(locs_GC_cropped,fullfile(path_output_ripley,condition,'GC',strcat(area_token,'_',condition,'_',num2str(synaptosomeID(j)),'.csv')),format)
            writeLocFile(locs_BC_cropped,fullfile(path_output_ripley,condition,'BC',strcat(area_token,'_',condition,'_',num2str(synaptosomeID(j)),'.csv')),format)

            X_RC = locs_RC_cropped.x;
            Y_RC = locs_RC_cropped.y;
            X_GC = locs_GC_cropped.x;
            Y_GC = locs_GC_cropped.y;
            X_BC = locs_BC_cropped.x;
            Y_BC = locs_BC_cropped.y;
                        
            num_Locs_RC_i(j,1) = size(X_RC,1);
            num_Locs_GC_i(j,1) = size(X_GC,1);
            num_Locs_BC_i(j,1) = size(X_BC,1);

            % Get images of cropped regions and write them away
            sigma = 150/magnification;
            intensities_RC = locs_RC_cropped.intensity;
            intensities_GC = locs_GC_cropped.intensity;
            intensities_BC = locs_BC_cropped.intensity;

            img_RC = generateImage(X_RC-min(X_RC)+1, Y_RC-min(Y_RC)+1, sigma, intensities_RC, windowsize/magnification, pixelsize, magnification);
            img_GC = generateImage(X_GC-min(X_GC)+1, Y_GC-min(Y_GC)+1, sigma, intensities_GC, windowsize/magnification, pixelsize, magnification);
            img_BC = generateImage(X_BC-min(X_BC)+1, Y_BC-min(Y_BC)+1, sigma, intensities_BC, windowsize/magnification, pixelsize, magnification);

            img_name_RC = strcat(char(sampleIDs(i)),'_RC_synaptosome_',num2str(synaptosomeID(j)),'.png');
            img_name_GC = strcat(char(sampleIDs(i)),'_GC_synaptosome_',num2str(synaptosomeID(j)),'.png');
            img_name_BC = strcat(char(sampleIDs(i)),'_BC_synaptosome_',num2str(synaptosomeID(j)),'.png');

            imwrite(flip(flip(img_RC,1),2),fullfile(path_output_ripley,condition,'images',img_name_RC));
            imwrite(flip(flip(img_GC,1),2),fullfile(path_output_ripley,condition,'images',img_name_GC));
            imwrite(flip(flip(img_BC,1),2),fullfile(path_output_ripley,condition,'images',img_name_BC));

            % Ripley's K single colour analysis
            [~,~,K_R,~,~] = CalcRipleyK(X_RC,Y_RC, X_RC,Y_RC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_G,~,~] = CalcRipleyK(X_GC,Y_GC, X_GC,Y_GC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_B,~,~] = CalcRipleyK(X_BC,Y_BC, X_BC,Y_BC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            ripley_RR_i(:,j) = sqrt(K_R'/pi) - r_hist;
            ripley_GG_i(:,j) = sqrt(K_G'/pi) - r_hist;
            ripley_BB_i(:,j) = sqrt(K_B'/pi) - r_hist;

            % Ripley's K two colour analysis
            [~,~,K_RG,~,~] = CalcRipleyK(X_RC,Y_RC, X_GC,Y_GC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_RB,~,~] = CalcRipleyK(X_RC,Y_RC, X_BC,Y_BC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            [~,~,K_GB,~,~] = CalcRipleyK(X_GC,Y_GC, X_BC,Y_BC, [Fov*pixelsize Fov*pixelsize], (windowsize*pixelsize)^2, R_max, r_step);
            ripley_RG_i(:,j) = sqrt(K_RG'/pi) - r_hist;
            ripley_RB_i(:,j) = sqrt(K_RB'/pi) - r_hist;
            ripley_GB_i(:,j) = sqrt(K_GB'/pi) - r_hist;
        end
        % Concatenate results ripley with previous samples
        ripley_RR = [ripley_RR ripley_RR_i];
        ripley_GG = [ripley_GG ripley_GG_i];
        ripley_BB = [ripley_BB ripley_BB_i];
        ripley_RG = [ripley_RG ripley_RG_i];
        ripley_RB = [ripley_RB ripley_RB_i];
        ripley_GB = [ripley_GB ripley_GB_i];
        
        numLocs_RC_EGTAK = [numLocs_RC_EGTAK num_Locs_RC_i'];
        numLocs_GC_EGTAK = [numLocs_GC_EGTAK num_Locs_GC_i'];
        numLocs_BC_EGTAK = [numLocs_BC_EGTAK num_Locs_BC_i'];
        
        rejectionRate_RC_EGTAK = [rejectionRate_RC_EGTAK rejectionRate_RC_EGTAK_i'];
        rejectionRate_GC_EGTAK = [rejectionRate_GC_EGTAK rejectionRate_GC_EGTAK_i'];
        rejectionRate_BC_EGTAK = [rejectionRate_BC_EGTAK rejectionRate_BC_EGTAK_i'];
        
        
    end

    H_all_egtak_RC = ripley_RR;
    H_all_egtak_GC = ripley_GG;
    H_all_egtak_BC = ripley_BB;
    H_all_egtak_RG = ripley_RG;
    H_all_egtak_RB = ripley_RB;
    H_all_egtak_GB = ripley_GB;

    results_EGTAK.numLocs_RC = numLocs_RC_EGTAK';
    results_EGTAK.numLocs_GC = numLocs_GC_EGTAK';
    results_EGTAK.numLocs_BC = numLocs_BC_EGTAK';
    
    results_EGTAK.rejectionRate_RC = rejectionRate_RC_EGTAK';
    results_EGTAK.rejectionRate_GC = rejectionRate_GC_EGTAK';
    results_EGTAK.rejectionRate_BC = rejectionRate_BC_EGTAK';
    
    % Write away Ripley's functions before removing NaN columns (to retain
    % information about from which synaptosome the curve came)
    path_ripley_egtak = fullfile(path_output_ripley,'H_egtak.mat');
    save(path_ripley_egtak,'ripley_RR','ripley_GG','ripley_BB',...
                          'ripley_RG','ripley_RB','ripley_GB');

    % Write away Ripley's functions
    path_H = fullfile(path_output_ripley,'H_all.mat');
    save(path_H,'H_all_phys_RC','H_all_phys_GC','H_all_phys_BC',...
                'H_all_phys_RG','H_all_phys_RB','H_all_phys_GB',...
                'H_all_egta_RC','H_all_egta_GC','H_all_egta_BC',...
                'H_all_egta_RG','H_all_egta_RB','H_all_egta_GB',...
                'H_all_egtak_RC','H_all_egtak_GC','H_all_egtak_BC',...
                'H_all_egtak_RG','H_all_egtak_RB','H_all_egtak_GB');

    % clear unfiltered localisation tables to free memory
    clear locs_RC_prefilter locs_GC_prefilter locs_BC_prefilter;
    
    
    %% Get clustersize from Ripley's L(r) - r function

    % Get clustersize for mCling from Ripley's curves
    [clustersize_phys_RC,  Mean_phys_RC,  SD_phys_RC]  = getStatsClustersize(r_hist,H_all_phys_RC);
    [clustersize_egta_RC,  Mean_egta_RC,  SD_egta_RC]  = getStatsClustersize(r_hist,H_all_egta_RC);
    [clustersize_egtak_RC, Mean_egtak_RC, SD_egtak_RC] = getStatsClustersize(r_hist,H_all_egtak_RC);

    % Get clustersize for a-sunclein from Ripley's curves
    [clustersize_phys_GC,  Mean_phys_GC,  SD_phys_GC]  = getStatsClustersize(r_hist,H_all_phys_GC);
    [clustersize_egta_GC,  Mean_egta_GC,  SD_egta_GC]  = getStatsClustersize(r_hist,H_all_egta_GC);
    [clustersize_egtak_GC, Mean_egtak_GC, SD_egtak_GC] = getStatsClustersize(r_hist,H_all_egtak_GC);

    % Get clustersize for VAMP2 from Ripley's curves
    [clustersize_phys_BC,  Mean_phys_BC,  SD_phys_BC]  = getStatsClustersize(r_hist,H_all_phys_BC);
    [clustersize_egta_BC,  Mean_egta_BC,  SD_egta_BC]  = getStatsClustersize(r_hist,H_all_egta_BC);
    [clustersize_egtak_BC, Mean_egtak_BC, SD_egtak_BC] = getStatsClustersize(r_hist,H_all_egtak_BC);

    % Add the values to the results tables
    results_PHYS.clustersizeRC = clustersize_phys_RC;
    results_PHYS.clustersizeGC = clustersize_phys_GC;
    results_PHYS.clustersizeBC = clustersize_phys_BC;

    results_EGTA.clustersizeRC = clustersize_egta_RC;
    results_EGTA.clustersizeGC = clustersize_egta_GC;
    results_EGTA.clustersizeBC = clustersize_egta_BC;

    results_EGTAK.clustersizeRC = clustersize_egtak_RC;
    results_EGTAK.clustersizeGC = clustersize_egtak_GC;
    results_EGTAK.clustersizeBC = clustersize_egtak_BC;

    path_clustersizes = fullfile(path_output_ripley,'clustersizes.mat');
    save(path_clustersizes,'clustersize_phys_RC','clustersize_egta_RC','clustersize_egtak_RC',...
                           'clustersize_phys_GC','clustersize_egta_GC','clustersize_egtak_GC',...
                           'clustersize_phys_BC','clustersize_egta_BC','clustersize_egtak_BC');

    % Convert results to tables and write away as tab-delimited txt file
    % (so they can be easily read in in R)

    % Red channel
    clear var phys egta egtak
    [phys{ 1:size(clustersize_phys_RC,1)}]  = deal('phys');
    [egta{ 1:size(clustersize_egta_RC,1)}]  = deal('egta');
    [egtak{1:size(clustersize_egtak_RC,1)}] = deal('egtak');
    results_clustersize_RC = array2table([clustersize_phys_RC;clustersize_egta_RC;clustersize_egtak_RC],'VariableNames',{'clustersize'});
    results_clustersize_RC.conditions = [phys egta egtak]';
    results_clustersize_RC = results_clustersize_RC(:,[end 1:end-1]);
    path_results_clustersize_RC = fullfile(path_output_ripley,'clustersizes_RC.txt');
    writetable(results_clustersize_RC,path_results_clustersize_RC,'Delimiter','\t');

    % Green channel
    clear var phys egta egtak
    [phys{ 1:size(clustersize_phys_GC,1)}]  = deal('phys');
    [egta{ 1:size(clustersize_egta_GC,1)}]  = deal('egta');
    [egtak{1:size(clustersize_egtak_GC,1)}] = deal('egtak');
    results_clustersize_GC = array2table([clustersize_phys_GC;clustersize_egta_GC;clustersize_egtak_GC],'VariableNames',{'clustersize'});
    results_clustersize_GC.conditions = [phys egta egtak]';
    results_clustersize_GC = results_clustersize_GC(:,[end 1:end-1]);
    path_results_clustersize_GC = fullfile(path_output_ripley,'clustersizes_GC.txt');
    writetable(results_clustersize_GC,path_results_clustersize_GC,'Delimiter','\t');

    % Blue channel
    clear var phys egta egtak
    [phys{ 1:size(clustersize_phys_BC,1)}]  = deal('phys');
    [egta{ 1:size(clustersize_egta_BC,1)}]  = deal('egta');
    [egtak{1:size(clustersize_egtak_BC,1)}] = deal('egtak');
    results_clustersize_BC = array2table([clustersize_phys_BC;clustersize_egta_BC;clustersize_egtak_BC],'VariableNames',{'clustersize'});
    results_clustersize_BC.conditions = [phys egta egtak]';
    results_clustersize_BC = results_clustersize_BC(:,[end 1:end-1]);
    path_results_clustersize_BC = fullfile(path_output_ripley,'clustersizes_BC.txt');
    writetable(results_clustersize_BC,path_results_clustersize_BC,'Delimiter','\t');


    %% Get intercluster distances from Ripley's L(r) - r function

    % Get clustersize for mCling from Ripley's curves
    [interclusterdist_phys_RG,  Mean_phys_RG,  SD_phys_RG]  = getStatsClustersize(r_hist,H_all_phys_RG);
    [interclusterdist_egta_RG,  Mean_egta_RG,  SD_egta_RG]  = getStatsClustersize(r_hist,H_all_egta_RG);
    [interclusterdist_egtak_RG, Mean_egtak_RG, SD_egtak_RG] = getStatsClustersize(r_hist,H_all_egtak_RG);

    % Get clustersize for a-sunclein from Ripley's curves
    [interclusterdist_phys_RB,  Mean_phys_RB,  SD_phys_RB]  = getStatsClustersize(r_hist,H_all_phys_RB);
    [interclusterdist_egta_RB,  Mean_egta_RB,  SD_egta_RB]  = getStatsClustersize(r_hist,H_all_egta_RB);
    [interclusterdist_egtak_RB, Mean_egtak_RB, SD_egtak_RB] = getStatsClustersize(r_hist,H_all_egtak_RB);

    % Get clustersize for VAMP2 from Ripley's curves
    [interclusterdist_phys_GB,  Mean_phys_GB,  SD_phys_GB]  = getStatsClustersize(r_hist,H_all_phys_GB);
    [interclusterdist_egta_GB,  Mean_egta_GB,  SD_egta_GB]  = getStatsClustersize(r_hist,H_all_egta_GB);
    [interclusterdist_egtak_GB, Mean_egtak_GB, SD_egtak_GB] = getStatsClustersize(r_hist,H_all_egtak_GB);


    % Add the values to the results tables
    results_PHYS.interclusterdistRG = interclusterdist_phys_RG;
    results_PHYS.interclusterdistRB = interclusterdist_phys_RB;
    results_PHYS.interclusterdistGB = interclusterdist_phys_GB;

    results_EGTA.interclusterdistRG = interclusterdist_egta_RG;
    results_EGTA.interclusterdistRB = interclusterdist_egta_RB;
    results_EGTA.interclusterdistGB = interclusterdist_egta_GB;

    results_EGTAK.interclusterdistRG = interclusterdist_egtak_RG;
    results_EGTAK.interclusterdistRB = interclusterdist_egtak_RB;
    results_EGTAK.interclusterdistGB = interclusterdist_egtak_GB;

    path_clustersizes = fullfile(path_output_ripley,'interclusterdist.mat');
    save(path_clustersizes,'interclusterdist_phys_RG','interclusterdist_egta_RG','interclusterdist_egtak_RG',...
                           'interclusterdist_phys_RB','interclusterdist_egta_RB','interclusterdist_egtak_RB',...
                           'interclusterdist_phys_GB','interclusterdist_egta_GB','interclusterdist_egtak_GB');

    % Convert results to tables and write away as tab-delimited txt file
    % (so they can be easily read in in R)

    % Red channel
    clear var phys egta egtak
    [phys{ 1:size(interclusterdist_phys_RG,1)}]  = deal('phys');
    [egta{ 1:size(interclusterdist_egta_RG,1)}]  = deal('egta');
    [egtak{1:size(interclusterdist_egtak_RG,1)}] = deal('egtak');
    results_interclusterdist_RG = array2table([interclusterdist_phys_RG;interclusterdist_egta_RG;interclusterdist_egtak_RG],'VariableNames',{'clustersize'});
    results_interclusterdist_RG.conditions = [phys egta egtak]';
    results_interclusterdist_RG = results_interclusterdist_RG(:,[end 1:end-1]);
    path_results_interclusterdist_RG = fullfile(path_output_ripley,'interclusterdist_RG.txt');
    writetable(results_interclusterdist_RG,path_results_interclusterdist_RG,'Delimiter','\t');

    % Green channel
    clear var phys egta egtak
    [phys{ 1:size(interclusterdist_phys_RB,1)}]  = deal('phys');
    [egta{ 1:size(interclusterdist_egta_RB,1)}]  = deal('egta');
    [egtak{1:size(interclusterdist_egtak_RB,1)}] = deal('egtak');
    results_interclusterdist_RB = array2table([interclusterdist_phys_RB;interclusterdist_egta_RB;interclusterdist_egtak_RB],'VariableNames',{'clustersize'});
    results_interclusterdist_RB.conditions = [phys egta egtak]';
    results_interclusterdist_RB = results_interclusterdist_RB(:,[end 1:end-1]);
    path_results_interclusterdist_RB = fullfile(path_output_ripley,'interclusterdist_RB.txt');
    writetable(results_interclusterdist_RB,path_results_interclusterdist_RB,'Delimiter','\t');

    % Blue channel
    clear var phys egta egtak
    [phys{ 1:size(interclusterdist_phys_GB,1)}]  = deal('phys');
    [egta{ 1:size(interclusterdist_egta_GB,1)}]  = deal('egta');
    [egtak{1:size(interclusterdist_egtak_GB,1)}] = deal('egtak');
    results_interclusterdist_GB = array2table([interclusterdist_phys_GB;interclusterdist_egta_GB;interclusterdist_egtak_GB],'VariableNames',{'clustersize'});
    results_interclusterdist_GB.conditions = [phys egta egtak]';
    results_interclusterdist_GB = results_interclusterdist_GB(:,[end 1:end-1]);
    path_results_interclusterdist_GB = fullfile(path_output_ripley,'interclusterdist_GB.txt');
    writetable(results_interclusterdist_GB,path_results_interclusterdist_GB,'Delimiter','\t');


    %% Merge tables of the three conditions after filtering

    results_combined_after_filtering = [results_PHYS; results_EGTA; results_EGTAK];

    clear var phys egta egtak
    [phys{1:size(results_PHYS,1)}] = deal('phys');
    [egta{1:size(results_EGTA,1)}] = deal('egta');
    [egtak{1:size(results_EGTAK,1)}] = deal('egtak');
    condition_column = [phys egta egtak];

    results_combined_after_filtering.condition = condition_column';
    results_combined_after_filtering = results_combined_after_filtering(:,[end 1:end-1]);

    % Remove rows that have -1 value in any of the Ripley's columns
    results_combined_after_filtering(ismember(results_combined_after_filtering.clustersizeRC,-1),:)=[];
    results_combined_after_filtering(ismember(results_combined_after_filtering.clustersizeGC,-1),:)=[];
    results_combined_after_filtering(ismember(results_combined_after_filtering.clustersizeBC,-1),:)=[];
    results_combined_after_filtering(ismember(results_combined_after_filtering.interclusterdistRG,-1),:)=[];
    results_combined_after_filtering(ismember(results_combined_after_filtering.interclusterdistRB,-1),:)=[];
    results_combined_after_filtering(ismember(results_combined_after_filtering.interclusterdistGB,-1),:)=[];

    
    %% Concatenate results to that of other repeats (if any)
    
    % Create an extra column for repeat
    clear var repeat
    [repeat{1:size(results_combined_after_filtering,1)}] = deal(repeats{n});
    results_combined_after_filtering.Repeat = repeat';
    results_combined_after_filtering = results_combined_after_filtering(:,[end 1:end-1]);
    
    if ~exist('results_repeats_pooled')
        results_repeats_pooled = results_combined_after_filtering;
    else
        results_repeats_pooled = [results_repeats_pooled; results_combined_after_filtering];
    end
    
    %% Save results

    path_results_combined = fullfile(path_output,'results_combined_after_overlap_threshold.mat');
    save(path_results_combined,'results_combined_after_filtering');

    path_results_combined = fullfile(path_output,'results_combined_after_overlap_threshold.txt');
    writetable(results_combined_after_filtering,path_results_combined,'Delimiter','\t');


    %% Mark detected synaptosomes on three-colour reconstructions

    path_detected = fullfile(path_output,'detected_synaptosomes');
    mkdir(path_detected);

    % PHYS condition
    results_PHYS = results_combined_after_filtering(...
        strcmp(results_combined_after_filtering.condition,'phys'),:);
    sampleIDs = unique(results_PHYS.sampleID);
    for i = 1:length(sampleIDs)

        % Read in image
        currentSample = sampleIDs(i);
        area_token = strsplit(char(currentSample),'_');
        area_token = area_token{1};
        filename = strcat(area_token,imgtype);
        img = imread(fullfile(dir_PHYS,filename));
        img = flip(flip(img,1),2);

        % Get coordinates of synaptosomes
        results_PHYS_i = results_PHYS(strcmp(results_PHYS.sampleID,currentSample),:);
        centroids = [results_PHYS_i.yCentroid results_PHYS_i.xCentroid];

        % Draw circles around coordinates
        img = drawCircles(img, centroids, radius, colour);
        if show; figure; imshow(img); end
        img_name = strcat(char(currentSample),'_detected.png');
        if flagprint; imwrite(img, fullfile(path_detected,img_name)); end
    end

    % EGTA condition
    results_EGTA = results_combined_after_filtering(...
        strcmp(results_combined_after_filtering.condition,'egta'),:);
    sampleIDs = unique(results_EGTA.sampleID);
    for i = 1:length(sampleIDs)

        % Read in image
        currentSample = sampleIDs(i);
        area_token = strsplit(char(currentSample),'_');
        area_token = area_token{1};
        filename = strcat(area_token,imgtype);
        img = imread(fullfile(dir_EGTA,filename));
        img = flip(flip(img,1),2);

        % Get coordinates of synaptosomes
        results_EGTA_i = results_EGTA(strcmp(results_EGTA.sampleID,currentSample),:);
        centroids = [results_EGTA_i.yCentroid results_EGTA_i.xCentroid];

        % Draw circles around coordinates
        img = drawCircles(img, centroids, radius, colour);
        if show; figure; imshow(img); end
        img_name = strcat(char(currentSample),'_detected.png');
        if flagprint; imwrite(img, fullfile(path_detected,img_name)); end
    end

    % EGTAK condition
    results_EGTAK = results_combined_after_filtering(...
        strcmp(results_combined_after_filtering.condition,'egtak'),:);
    sampleIDs = unique(results_EGTAK.sampleID);
    for i = 1:length(sampleIDs)

        % Read in image
        currentSample = sampleIDs(i);
        area_token = strsplit(char(currentSample),'_');
        area_token = area_token{1};
        filename = strcat(area_token,imgtype);
        img = imread(fullfile(dir_EGTAK,filename));
        img = flip(flip(img,1),2);

        % Get coordinates of synaptosomes
        results_EGTAK_i = results_EGTAK(strcmp(results_EGTAK.sampleID,currentSample),:);
        centroids = [results_EGTAK_i.yCentroid results_EGTAK_i.xCentroid];

        % Draw circles around coordinates
        img = drawCircles(img, centroids, radius, colour);
        if show; figure; imshow(img); end
        img_name = strcat(char(currentSample),'_detected.png');
        if flagprint; imwrite(img, fullfile(path_detected,img_name)); end
    end
    
    % Plot results of current repeat
    if plot_results_individual_repeats
        plotResults(path_output,results_PHYS,results_EGTA,results_EGTAK,magnification,flagprint)
    end
end

%% Write away pooled results from all the repeats in a new folder

% Create new output folder
path_output = fullfile(output_dir,'Results_combined','Repeats_combined');
mkdir(path_output);

% Write away results
path_results_combined = fullfile(path_output,'results_combined_after_overlap_threshold.mat');
save(path_results_combined,'results_repeats_pooled');
path_results_combined = fullfile(path_output,'results_combined_after_overlap_threshold.txt');
writetable(results_repeats_pooled,path_results_combined,'Delimiter','\t');

% Split results table into conditions (for convenience when plotting)
results_PHYS  = results_repeats_pooled(strcmp(results_repeats_pooled.condition,'phys'),:);
results_EGTA  = results_repeats_pooled(strcmp(results_repeats_pooled.condition,'egta'),:);
results_EGTAK = results_repeats_pooled(strcmp(results_repeats_pooled.condition,'egtak'),:);

% Plot all results
plotResults(path_output,results_PHYS,results_EGTA,results_EGTAK,magnification,flagprint)

diary off
toc
