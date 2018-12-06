clear all
close all
clc

% dir_PHYS  = 'F:\synaptosomes\Results_phys';
% dir_EGTA  = 'F:\synaptosomes\Results_egta';
% dir_EGTAK = 'F:\synaptosomes\Results_egtak';

 dir_PHYS  = 'E:\Experiments\synaptosomes\Results synaptosome_2nd_round\Results_phys';
 dir_EGTA  = 'E:\Experiments\synaptosomes\Results synaptosome_2nd_round\Results_egta';
 dir_EGTAK = 'E:\Experiments\synaptosomes\Results synaptosome_2nd_round\Results_egtak';

<<<<<<< HEAD

%output_dir = 'F:\Data\Synaptosomes\Experiment_4C\Results';
%output_dir = '/Volumes/WD Ezra/Dump';
%output_dir = fullfile(pwd,'testdata');
output_dir = fullfile('E:\Experiments\synaptosomes\Results synaptosome_2nd_round',filesep);
=======
% dir_PHYS  = '/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_phys';
% dir_EGTA  = '/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_egta';
% dir_EGTAK = '/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_egtak';

% dir_PHYS  = '/Volumes/WD Ezra/Data/Synaptosomes/Experiment_4C/Results/Results_phys';
% dir_EGTA  = '/Volumes/WD Ezra/Data/Synaptosomes/Experiment_4C/Results/Results_egta';
% dir_EGTAK = '/Volumes/WD Ezra/Data/Synaptosomes/Experiment_4C/Results/Results_egtak';

dir_PHYS  = fullfile(pwd,'testdata/Results_phys');
dir_EGTA  = fullfile(pwd,'testdata/Results_egta');
dir_EGTAK = fullfile(pwd,'testdata/Results_egtak');

%output_dir = 'F:\Data\Synaptosomes\Experiment_4C\Results';
%output_dir = '/Volumes/WD Ezra/Dump';
output_dir = fullfile(pwd,'testdata');
>>>>>>> upstream/master

format = 'thunderstorm';
magnification = 10;

imgtype = '_reconstruction_RGB.tif'; % image to mark detected synaptosomes on with circles
radius = 100; % radius of the circles
colour = [100 100 100]; % colour of the circles (triplet)

% Parameters for filtering on overlap between channels to detect synaptosomes
<<<<<<< HEAD
filterOverlapGreen = 1;  % 1 to filter on overlap between red and green to detect synaptosomes
filterOverlapBlue  = 1;  % 1 to filter on overlap between red and blue  to detect synaptosomes
minOverlapGreen    = 10;  % minimum overlap between red and green to detect synaptosomes
minOverlapBlue     = 10; % minimum overlap between red and green to detect synaptosomes

% Minimum distance between two synaptosomes
min_dist_between_synaptosomes = 300; % in nm
=======
filterOverlapGreen = 0;  % 1 to filter on overlap between red and green to detect synaptosomes
filterOverlapBlue  = 1;  % 1 to filter on overlap between red and blue  to detect synaptosomes
minOverlapGreen    = 0;  % minimum overlap between red and green to detect synaptosomes
minOverlapBlue     = 20; % minimum overlap between red and green to detect synaptosomes

% Minimum distance between two synaptosomes
min_dist_between_synaptosomes = 1500; % in nm
>>>>>>> upstream/master

r_step = 10;
R_max = 1000;
windowsize = 80; % pixels
pixelsize = 117; % nm

show = 0; % to show intermediate results
flagprint  = 1; % set to 1 to save the fig visualization of the results


%% Creating output folder in output_dir

% Create new output folder
path_output = fullfile(output_dir,'Results_combined');
if exist(path_output, 'dir')
    opts.Interpreter = 'tex';
    opts.Default = 'Continue';
    quest = '\fontsize{12}An output folder ''Results combined'' already exists. If you continue, data in this folder might be overwritten.';
    answer = questdlg(quest,'Message','Cancel','Continue',opts);
    if strcmp(answer,'Continue')
        mkdir(path_output);
    else
        return
    end
else
    mkdir(path_output);
end

save(fullfile(path_output,'parameters.mat'));
diary(fullfile(path_output,'command_window_output.txt'));
 
%% Read in results of different ROIs and merge them

% PHYS --------------------------------------------------------------------

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
results_PHYS = array2table(results_PHYS,'VariableNames',{'synaptosomeID','xCentroid','yCentroid','Area','OverlapWithGreen','OverlapWithBlue','WeightedOverlapWithGreen','WeightedOverlapWithBlue'});

% Add overlapping green and blue area
results_PHYS.AreaGreen = (results_PHYS.OverlapWithGreen.*results_PHYS.Area)/100;
results_PHYS.AreaBlue  = (results_PHYS.OverlapWithBlue.*results_PHYS.Area)/100;

% Add sampleID column in front
results_PHYS.sampleID  = sampleIDs';
results_PHYS = results_PHYS(:,[end 1:end-1]);


% EGTA --------------------------------------------------------------------

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
results_EGTA = array2table(results_EGTA,'VariableNames',{'synaptosomeID','xCentroid','yCentroid','Area','OverlapWithGreen','OverlapWithBlue','WeightedOverlapWithGreen','WeightedOverlapWithBlue'});

% Add overlapping green and blue area
results_EGTA.AreaGreen = (results_EGTA.OverlapWithGreen.*results_EGTA.Area)/100;
results_EGTA.AreaBlue  = (results_EGTA.OverlapWithBlue.*results_EGTA.Area)/100;

% Add sampleID column in front
results_EGTA.sampleID  = sampleIDs';
results_EGTA = results_EGTA(:,[end 1:end-1]);


% EGTAK -------------------------------------------------------------------

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
results_EGTAK = array2table(results_EGTAK,'VariableNames',{'synaptosomeID','xCentroid','yCentroid','Area','OverlapWithGreen','OverlapWithBlue','WeightedOverlapWithGreen','WeightedOverlapWithBlue'});

% Add overlapping green and blue area
results_EGTAK.AreaGreen = (results_EGTAK.OverlapWithGreen.*results_EGTAK.Area)/100;
results_EGTAK.AreaBlue  = (results_EGTAK.OverlapWithBlue.*results_EGTAK.Area)/100;

% Add sampleID column in front
results_EGTAK.sampleID  = sampleIDs';
results_EGTAK = results_EGTAK(:,[end 1:end-1]);


%% Merge tables of the three conditions before filtering

results_combined = [results_PHYS; results_EGTA; results_EGTAK];

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


%% Detect synaptosomes by filtering on overlap (remove false positives )
% Many of the blobs in the red channel for which overlaps are calculated
% are not actually synaptosomes. If they don'toverlap with blobs in the
% blue channel, they probably aren't synaptosomes. So, we filter out all
% results for red blobs that don't show significant overlap with the green
% channel.

<<<<<<< HEAD
synapto_number_prefilter_PHYS = size(results_PHYS);
synapto_number_prefilter_EGTA = size(results_EGTA);
synapto_number_prefilter_EGTAK = size(results_EGTAK);




=======
>>>>>>> upstream/master
if filterOverlapGreen
    
    results_PHYS  = results_PHYS(results_PHYS.OverlapWithGreen   > minOverlapGreen,:);
    results_EGTA  = results_EGTA(results_EGTA.OverlapWithGreen   > minOverlapGreen,:);
    results_EGTAK = results_EGTAK(results_EGTAK.OverlapWithGreen > minOverlapGreen,:);
<<<<<<< HEAD

elseif filterOverlapBlue

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


%% Remove synaptosomes that are too close together

% Make a copy of the results before the proximity filtering
% (so only the overlap-filtered results)
results_PHYS_before_proximity_filter = results_PHYS;
results_EGTA_before_proximity_filter = results_EGTA;
results_EGTAK_before_proximity_filter = results_EGTAK;


% PHYS condition ----------------------------------------------------------

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

% -------------------------------------------------------------------------


% EGTA condition ----------------------------------------------------------

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

% -------------------------------------------------------------------------

% EGTA condition ----------------------------------------------------------

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

disp(['Synaptosomes PHYS post proximity filter: ' num2str(synapto_number_post_proximityFilter_PHYS(1))]);
disp(['Synaptosomes EGTA post proximity filter: ' num2str(synapto_number_post_proximityFilter_EGTA(1))]);
disp(['Synaptosomes EGTAK post proximity filter: ' num2str(synapto_number_post_proximityFilter_EGTAK(1))]);
=======

elseif filterOverlapBlue
>>>>>>> upstream/master

    results_PHYS  = results_PHYS(results_PHYS.OverlapWithBlue   > minOverlapBlue,:);
    results_EGTA  = results_EGTA(results_EGTA.OverlapWithBlue   > minOverlapBlue,:);
    results_EGTAK = results_EGTAK(results_EGTAK.OverlapWithBlue > minOverlapBlue,:);

end

%% Ripley's K analysis

% Create new subfolder in output folder
path_output_ripley = fullfile(path_output,'ripley');
mkdir(path_output_ripley);
mkdir(fullfile(path_output_ripley,'phys'));
mkdir(fullfile(path_output_ripley,'egta'));
mkdir(fullfile(path_output_ripley,'egtak'));


% PHYS condition ----------------------------------------------------------

condition = 'phys';

% Initialize
ripley_RR = [];
ripley_GG = [];
ripley_BB = [];
ripley_RG = [];
ripley_RB = [];
ripley_GB = [];

% Create subfolders
mkdir(fullfile(path_output_ripley,condition,'RC'));
mkdir(fullfile(path_output_ripley,condition,'GC'));
mkdir(fullfile(path_output_ripley,condition,'BC'));
mkdir(fullfile(path_output_ripley,condition,'images'));

% Get sample IDs and loop over them
sampleIDs = unique(results_PHYS.sampleID);
for i = 1:length(sampleIDs)
   
    % Load localisation files of sample
    area_token = strsplit(char(sampleIDs(i)),'_');
    area_token = area_token{1};
    locs_RC = readLocFile(fullfile(dir_PHYS,strcat(area_token,'_locs_RC_filtered.csv')),format);
    locs_GC = readLocFile(fullfile(dir_PHYS,strcat(area_token,'_locs_GC_filtered.csv')),format);
    locs_BC = readLocFile(fullfile(dir_PHYS,strcat(area_token,'_locs_BC_filtered.csv')),format);
    
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
    
    % Loop over all synaptosomes
    for j = 1:num_synaptosomes
        
        % Get coordinates of synaptosome
        x_c = X_c(j)*magnification;
        y_c = Y_c(j)*magnification;

        % Crop locfiles to region around synaptosome
        locs_RC_cropped = cropLocsAroundCentroid(locs_RC,x_c,y_c,(windowsize/magnification)*pixelsize);
        locs_GC_cropped = cropLocsAroundCentroid(locs_GC,x_c,y_c,(windowsize/magnification)*pixelsize);
        locs_BC_cropped = cropLocsAroundCentroid(locs_BC,x_c,y_c,(windowsize/magnification)*pixelsize);
        
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
end


% H_all_phys_RC = ripley_RR(:,all(~isnan(ripley_RR)));
% H_all_phys_GC = ripley_GG(:,all(~isnan(ripley_GG)));
% H_all_phys_BC = ripley_BB(:,all(~isnan(ripley_BB)));
% H_all_phys_RG = ripley_RG(:,all(~isnan(ripley_RG)));
% H_all_phys_RB = ripley_RB(:,all(~isnan(ripley_RB)));
% H_all_phys_GB = ripley_GB(:,all(~isnan(ripley_GB)));

H_all_phys_RC = ripley_RR;
H_all_phys_GC = ripley_GG;
H_all_phys_BC = ripley_BB;
H_all_phys_RG = ripley_RG;
H_all_phys_RB = ripley_RB;
H_all_phys_GB = ripley_GB;


% Write away Ripley's functions before removing NaN columns (to retain
% information about from which synaptosome the curve came)
path_ripley_phys = fullfile(path_output_ripley,'H_phys.mat');
save(path_ripley_phys,'ripley_RR','ripley_GG','ripley_BB',...
                      'ripley_RG','ripley_RB','ripley_GB');


% EGTA condition ----------------------------------------------------------

condition = 'egta';

% Initialize
ripley_RR = [];
ripley_GG = [];
ripley_BB = [];
ripley_RG = [];
ripley_RB = [];
ripley_GB = [];

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
    
    % Loop over all synaptosomes
    for j = 1:num_synaptosomes
        
        % Get coordinates of synaptosome
        x_c = X_c(j)*magnification;
        y_c = Y_c(j)*magnification;

        % Crop locfiles to region around synaptosome
        locs_RC_cropped = cropLocsAroundCentroid(locs_RC,x_c,y_c,(windowsize/magnification)*pixelsize);
        locs_GC_cropped = cropLocsAroundCentroid(locs_GC,x_c,y_c,(windowsize/magnification)*pixelsize);
        locs_BC_cropped = cropLocsAroundCentroid(locs_BC,x_c,y_c,(windowsize/magnification)*pixelsize);
        
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
end

% H_all_egta_RC = ripley_RR(:,all(~isnan(ripley_RR)));
% H_all_egta_GC = ripley_GG(:,all(~isnan(ripley_GG)));
% H_all_egta_BC = ripley_BB(:,all(~isnan(ripley_BB)));
% H_all_egta_RG = ripley_RG(:,all(~isnan(ripley_RG)));
% H_all_egta_RB = ripley_RB(:,all(~isnan(ripley_RB)));
% H_all_egta_GB = ripley_GB(:,all(~isnan(ripley_GB)));

H_all_egta_RC = ripley_RR;
H_all_egta_GC = ripley_GG;
H_all_egta_BC = ripley_BB;
H_all_egta_RG = ripley_RG;
H_all_egta_RB = ripley_RB;
H_all_egta_GB = ripley_GB;


% Write away Ripley's functions before removing NaN columns (to retain
% information about from which synaptosome the curve came)
path_ripley_egta = fullfile(path_output_ripley,'H_egta.mat');
save(path_ripley_egta,'ripley_RR','ripley_GG','ripley_BB',...
                      'ripley_RG','ripley_RB','ripley_GB');
                  

% EGTAK condition ---------------------------------------------------------

condition = 'egtak';

% Initialize
ripley_RR = [];
ripley_GG = [];
ripley_BB = [];
ripley_RG = [];
ripley_RB = [];
ripley_GB = [];

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
    
    % Loop over all synaptosomes
    for j = 1:num_synaptosomes
        
        % Get coordinates of synaptosome
        x_c = X_c(j)*magnification;
        y_c = Y_c(j)*magnification;

        % Crop locfiles to region around synaptosome
        locs_RC_cropped = cropLocsAroundCentroid(locs_RC,x_c,y_c,(windowsize/magnification)*pixelsize);
        locs_GC_cropped = cropLocsAroundCentroid(locs_GC,x_c,y_c,(windowsize/magnification)*pixelsize);
        locs_BC_cropped = cropLocsAroundCentroid(locs_BC,x_c,y_c,(windowsize/magnification)*pixelsize);
        
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
end

% H_all_egtak_RC = ripley_RR(:,all(~isnan(ripley_RR)));
% H_all_egtak_GC = ripley_GG(:,all(~isnan(ripley_GG)));
% H_all_egtak_BC = ripley_BB(:,all(~isnan(ripley_BB)));
% H_all_egtak_RG = ripley_RG(:,all(~isnan(ripley_RG)));
% H_all_egtak_RB = ripley_RB(:,all(~isnan(ripley_RB)));
% H_all_egtak_GB = ripley_GB(:,all(~isnan(ripley_GB)));

H_all_egtak_RC = ripley_RR;
H_all_egtak_GC = ripley_GG;
H_all_egtak_BC = ripley_BB;
H_all_egtak_RG = ripley_RG;
H_all_egtak_RB = ripley_RB;
H_all_egtak_GB = ripley_GB;


% Write away Ripley's functions before removing NaN columns (to retain
% information about from which synaptosome the curve came)
path_ripley_egtak = fullfile(path_output_ripley,'H_egtak.mat');
save(path_ripley_egtak,'ripley_RR','ripley_GG','ripley_BB',...
                      'ripley_RG','ripley_RB','ripley_GB');

% Write away Ripley's functions after removing NaN columns
path_H = fullfile(path_output_ripley,'H_all.mat');
save(path_H,'H_all_phys_RC','H_all_phys_GC','H_all_phys_BC',...
            'H_all_phys_RG','H_all_phys_RB','H_all_phys_GB',...
            'H_all_egta_RC','H_all_egta_GC','H_all_egta_BC',...
            'H_all_egta_RG','H_all_egta_RB','H_all_egta_GB',...
            'H_all_egtak_RC','H_all_egtak_GC','H_all_egtak_BC',...
            'H_all_egtak_RG','H_all_egtak_RB','H_all_egtak_GB');


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


<<<<<<< HEAD
=======
%% Remove synaptosomes that are too close together

unique_conditions = unique(results_combined_after_filtering.condition);

% Loop over different conditions in the results
for i=1:size(unique_conditions)
    
    % Get only results of current condition
    condition_i = unique_conditions{i};
    results_condition_i = results_combined_after_filtering(...
        strcmp(results_combined_after_filtering.condition,condition_i),:);
    
    % Loop over sample_ID in the results
    sample_IDs = unique(results_condition_i.sampleID);
    for j=1:size(sample_IDs,1)
        sample_ID_j = sample_IDs{j};
        results_sampleID = results_condition_i(...
            strcmp(results_condition_i.sampleID,sample_ID_j),:);
        
        % Get indeces of synaptosomes in this sample that are too close
        % together
        points = [results_sampleID.xCentroid results_sampleID.yCentroid];
        pairwise_dist_matrix = pdist2(points,points);
        to_remove = (pairwise_dist_matrix < min_dist_between_synaptosomes) - eye(size(points,1));
        index_to_remove = sum(to_remove) > 0;
        disp(['To remove: ' num2str(sum(index_to_remove))]);

        % Remove the synaptosomes that are too close together
        filtered_results_sampleID = results_sampleID(index_to_remove == 0,:);
        if ~exist('final_results')
            final_results = filtered_results_sampleID;
        else
            final_results = [final_results; filtered_results_sampleID];
        end
    end
end

results_combined_after_filtering = final_results;

>>>>>>> upstream/master
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


%% Plot and compare results from overlap analysis

path_figures = fullfile(path_output,'figures');
mkdir(path_figures);
path_figures_overlap = fullfile(path_figures,'overlap_analysis');
mkdir(path_figures_overlap);

% Groups
[condition_PHYS{1:size(results_PHYS,1)}]   = deal('PHYS');
[condition_EGTA{1:size(results_EGTA,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(results_EGTAK,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];


% Area of synaptosomes ----------------------------------------------------

areaRed_PHYS = results_PHYS.Area;
areaRed_EGTA = results_EGTA.Area;
areaRed_EGTAK = results_EGTAK.Area;
areaRed = [areaRed_PHYS; areaRed_EGTA; areaRed_EGTAK];

% Beeswarm and boxplot
fig1 = figure;
subplot(121)
h = plotSpread({areaRed_PHYS*magnification^2,areaRed_EGTA*magnification^2,areaRed_EGTAK*magnification^2},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Area (nm^2)');
title('Synaptosome area');

subplot(122)
boxplot(areaRed*magnification^2,conditions);
ylabel('Area (nm^2)');
title('Synaptosome area');

if flagprint
    savefig(fig1,fullfile(path_figures_overlap,'area_synaptosomes.fig'))
end


% Area of a-synuclein -----------------------------------------------------

areaGreen_PHYS = results_PHYS.AreaGreen;
areaGreen_EGTA = results_EGTA.AreaGreen;
areaGreen_EGTAK = results_EGTAK.AreaGreen;
areaGreen = [areaGreen_PHYS; areaGreen_EGTA; areaGreen_EGTAK];

% Beeswarm and boxplot
fig4 = figure;
subplot(121)
h = plotSpread({areaGreen_PHYS*magnification^2,areaGreen_EGTA*magnification^2,areaGreen_EGTAK*magnification^2},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Area (nm^2)');
title('a-synuclein area');

subplot(122)
boxplot(areaGreen*magnification^2,conditions);
ylabel('Area (nm^2)');
title('a-synuclein area');

if flagprint
    savefig(fig4,fullfile(path_figures_overlap,'area-a-synuclein.fig'))
end


% Area of VAMP2 -----------------------------------------------------------

areaBlue_PHYS = results_PHYS.AreaBlue;
areaBlue_EGTA = results_EGTA.AreaBlue;
areaBlue_EGTAK = results_EGTAK.AreaBlue;
areaBlue = [areaBlue_PHYS; areaBlue_EGTA; areaBlue_EGTAK];

% Beeswarm and boxplot
fig5 = figure;
subplot(121)
h = plotSpread({areaBlue_PHYS*magnification^2,areaBlue_EGTA*magnification^2,areaBlue_EGTAK*magnification^2},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Area (nm^2)');
title('VAMP2 area');

subplot(122)
boxplot(areaBlue*magnification^2,conditions);
ylabel('Area (nm^2)');
title('VAMP2 area');

if flagprint
    savefig(fig5,fullfile(path_figures_overlap,'area-VAMP2.fig'))
end


% Overlap of synaptosomes with a-synuclein --------------------------------

overlapGreen_PHYS = results_PHYS.OverlapWithGreen;
overlapGreen_EGTA = results_EGTA.OverlapWithGreen;
overlapGreen_EGTAK = results_EGTAK.OverlapWithGreen;
overlapGreen = [overlapGreen_PHYS; overlapGreen_EGTA; overlapGreen_EGTAK];

% Beeswarm and boxplot
fig2 = figure;
subplot(121)
h = plotSpread({overlapGreen_PHYS,overlapGreen_EGTA,overlapGreen_EGTAK},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('% area overlap'); ylim([0 100]);
title({'Unweighted overlap','mCling and a-synuclein'});

subplot(122)
boxplot(overlapGreen,conditions);
ylabel('% area overlap'); ylim([0 100]);
title({'Unweighted overlap','mCling and a-synuclein'});

if flagprint
    savefig(fig2,fullfile(path_figures_overlap,'unweighted-overlap-mcling-a-synuclein.fig'))
    % export figure to a pdf for display
    filename = fullfile(path_figures_overlap,'unweighted-overlap-mcling-a-synuclein');
    print(fig2,filename,'-dpdf','-r300')
end


% Overlap of synaptosomes with VAMP2 --------------------------------------

overlapBlue_PHYS = results_PHYS.OverlapWithBlue;
overlapBlue_EGTA = results_EGTA.OverlapWithBlue;
overlapBlue_EGTAK = results_EGTAK.OverlapWithBlue;
overlapBlue = [overlapBlue_PHYS; overlapBlue_EGTA; overlapBlue_EGTAK];

% Beeswarm and boxplot
fig3 = figure;
subplot(121)
h = plotSpread({overlapBlue_PHYS,overlapBlue_EGTA,overlapBlue_EGTAK},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('% area overlap'); ylim([0 100]);
title({'Unweighted overlap','mCling and VAMP2'});

subplot(122)
boxplot(overlapBlue,conditions);
ylabel('% area overlap'); ylim([0 100]);
title({'Unweighted overlap','mCling and VAMP2'});

if flagprint
    savefig(fig3,fullfile(path_figures_overlap,'unweighted-overlap-mcling-vamp2.fig'))
    filename = fullfile(path_figures_overlap,'unweighted-overlap-mcling-vamp2');
    print(fig3,filename,'-dpdf','-r300')
end


% Weighted overlap of synaptosomes with a-synuclein -----------------------

weightedOverlapWithGreen_PHYS  = results_PHYS.WeightedOverlapWithGreen;
weightedOverlapWithGreen_EGTA  = results_EGTA.WeightedOverlapWithGreen;
weightedOverlapWithGreen_EGTAK = results_EGTAK.WeightedOverlapWithGreen;
weightedOverlapWithGreen = [weightedOverlapWithGreen_PHYS;
                            weightedOverlapWithGreen_EGTA;
                            weightedOverlapWithGreen_EGTAK];
fig6 = figure;
subplot(121)
h = plotSpread({weightedOverlapWithGreen_PHYS,weightedOverlapWithGreen_EGTA,weightedOverlapWithGreen_EGTAK},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Weighted overlap');
title({'Weighted overlap','mCling and a-synuclein'});

subplot(122)
boxplot(weightedOverlapWithGreen,conditions);
ylabel('Weighted overlap');
title({'Weighted overlap','mCling and a-synuclein'});

if flagprint
    savefig(fig6,fullfile(path_figures_overlap,'weighted-overlap-mcling-a-synuclein.fig'))
    filename = fullfile(path_figures_overlap,'weighted-overlap-mcling-a-synuclein');
    print(fig6,filename,'-dpdf','-r300')
end


% Weighted overlap of synaptosomes with VAMP2 -----------------------------

weightedOverlapWithBlue_PHYS  = results_PHYS.WeightedOverlapWithBlue;
weightedOverlapWithBlue_EGTA  = results_EGTA.WeightedOverlapWithBlue;
weightedOverlapWithBlue_EGTAK = results_EGTAK.WeightedOverlapWithBlue;
weightedOverlapWithBlue = [weightedOverlapWithBlue_PHYS;
                           weightedOverlapWithBlue_EGTA;
                           weightedOverlapWithBlue_EGTAK];
fig7 = figure;
subplot(121)
h = plotSpread({weightedOverlapWithBlue_PHYS,weightedOverlapWithBlue_EGTA,weightedOverlapWithBlue_EGTAK},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Weighted overlap');
title({'Weighted overlap','mCling and VAMP2'});

subplot(122)
boxplot(weightedOverlapWithBlue,conditions);
ylabel('Weighted overlap');
title({'Weighted overlap','mCling and VAMP2'});

if flagprint
    savefig(fig7,fullfile(path_figures_overlap,'weighted-overlap-mcling-vamp2.fig'))
    filename = fullfile(path_figures_overlap,'weighted-overlap-mcling-vamp2');
    print(fig7,filename,'-dpdf','-r300')
end

%% Plot and compare results from Ripley analysis

path_figures_ripley = fullfile(path_figures,'ripley_analysis');
mkdir(path_figures_ripley);

% Clustersize mCling from Ripley ------------------------------------------

% Get groups
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(clustersize_phys_RC,1)}]   = deal('PHYS');
[condition_EGTA{1:size(clustersize_egta_RC,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(clustersize_egtak_RC,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot clustersizes estimated from Ripley's curves
fig8 = figure;

subplot(121) % beeswarm plot
h = plotSpread({clustersize_phys_RC,clustersize_egta_RC,clustersize_egtak_RC},'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Cluster size (nm)');
title('Clustersize mCling');
set(gca,'fontsize',14);
%ylim([0 1000])

subplot(122) % boxplot
boxplot([clustersize_phys_RC; clustersize_egta_RC; clustersize_egtak_RC],conditions);
ylabel('Cluster size (nm)');
title('Clustersize mCling');
set(gca,'fontsize',14);
%ylim([0 1000])

if flagprint
    savefig(fig8,fullfile(path_figures_ripley,'clustersize_mCling.fig'))
end

% One-way ANOVA
p = anova1([clustersize_phys_RC; clustersize_egta_RC; clustersize_egtak_RC],conditions,'off');
if p < 0.05
    disp(['Significant difference in mCling clustersize between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in mCling clustersize between conditions (p-value = ' num2str(p) ')'])
end

% Clustersize a-synuclein from Ripley -------------------------------------

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(clustersize_phys_GC,1)}]   = deal('PHYS');
[condition_EGTA{1:size(clustersize_egta_GC,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(clustersize_egtak_GC,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot clustersizes estimated from Ripley's curves
fig9 = figure;
subplot(121) % beeswarm plot
h = plotSpread({clustersize_phys_GC,clustersize_egta_GC,clustersize_egtak_GC},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Cluster size (nm)');
title('Clustersize a-synuclein');
set(gca,'fontsize',14);
%ylim([0 1000])

subplot(122) % boxplot
boxplot([clustersize_phys_GC; clustersize_egta_GC; clustersize_egtak_GC],conditions);
ylabel('Cluster size (nm)');
title('Clustersize a-synuclein');
set(gca,'fontsize',14);
%ylim([0 1000])

if flagprint
    savefig(fig9,fullfile(path_figures_ripley,'clustersize_a-synuclein.fig'))
end

% One-way ANOVA
p = anova1([clustersize_phys_GC; clustersize_egta_GC; clustersize_egtak_GC],conditions,'off');
if p < 0.05
    disp(['Significant difference in a-synuclein clustersize between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in a-synuclein clustersize between conditions (p-value = ' num2str(p) ')'])
end

% Clustersize VAMP2 from Ripley -------------------------------------------

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(clustersize_phys_BC,1)}]   = deal('PHYS');
[condition_EGTA{1:size(clustersize_egta_BC,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(clustersize_egtak_BC,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot clustersizes estimated from Ripley's curves
fig10 = figure;
subplot(121) % beeswarm plot
h = plotSpread({clustersize_phys_BC,clustersize_egta_BC,clustersize_egtak_BC},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Cluster size (nm)');
title('Clustersize VAMP2');
set(gca,'fontsize',14);
%ylim([0 1000])

subplot(122) % boxplot
boxplot([clustersize_phys_BC; clustersize_egta_BC; clustersize_egtak_BC],conditions);
ylabel('Cluster size (nm)');
title('Clustersize VAMP2');
set(gca,'fontsize',14);
%ylim([0 1000])

if flagprint
    savefig(fig10,fullfile(path_figures_ripley,'clustersize_VAMP2.fig'))
end

% One-way ANOVA
p = anova1([clustersize_phys_BC; clustersize_egta_BC; clustersize_egtak_BC],conditions,'off');
if p < 0.05
    disp(['Significant difference in VAMP2 clustersize between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in VAMP2 clustersize between conditions (p-value = ' num2str(p) ')'])
end


% Intercluster distance between mCling and a-synuclein from Ripley --------

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(interclusterdist_phys_RG,1)}]   = deal('PHYS');
[condition_EGTA{1:size(interclusterdist_egta_RG,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(interclusterdist_egtak_RG,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot intercluster distances estimated from Ripley's curves
fig11 = figure;
subplot(121) % beeswarm plot
h = plotSpread({interclusterdist_phys_RG,interclusterdist_egta_RG,interclusterdist_egtak_RG},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','mCling and a-synuclein'});
set(gca,'fontsize',14);

subplot(122) % boxplot
boxplot([interclusterdist_phys_RG; interclusterdist_egta_RG; interclusterdist_egtak_RG],conditions);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','mCling and a-synuclein'});
set(gca,'fontsize',14);

if flagprint
    savefig(fig11,fullfile(path_figures_ripley,'intercluster_distance_mCling_a-synuclein.fig'))
end

% One-way ANOVA
p = anova1([interclusterdist_phys_RG; interclusterdist_egta_RG; interclusterdist_egtak_RG],conditions,'off');
if p < 0.05
    disp(['Significant difference in mCling/a-synuclein intercluster distance between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in mCling/a-synuclein intercluster distance between conditions (p-value = ' num2str(p) ')'])
end


% Intercluster distance between mCling and VAMP2 from Ripley --------------

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(interclusterdist_phys_RB,1)}]   = deal('PHYS');
[condition_EGTA{1:size(interclusterdist_egta_RB,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(interclusterdist_egtak_RB,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot intercluster distances estimated from Ripley's curves
fig12 = figure;
subplot(121) % beeswarm plot
h = plotSpread({interclusterdist_phys_RB,interclusterdist_egta_RB,interclusterdist_egtak_RB},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','mCling and VAMP2'});
set(gca,'fontsize',14);

subplot(122) % boxplot
boxplot([interclusterdist_phys_RB; interclusterdist_egta_RB; interclusterdist_egtak_RB],conditions);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','mCling and VAMP2'});
set(gca,'fontsize',14);

if flagprint
    savefig(fig12,fullfile(path_figures_ripley,'intercluster_distance_mCling_VAMP2.fig'))
end

% One-way ANOVA
p = anova1([interclusterdist_phys_RB; interclusterdist_egta_RB; interclusterdist_egtak_RB],conditions,'off');
if p < 0.05
    disp(['Significant difference in mCling/VAMP2 intercluster distance between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in mCling/VAMP2 intercluster distance between conditions (p-value = ' num2str(p) ')'])
end


% Intercluster distance between a-synuclein and VAMP2 from Ripley ---------

% Get groups (needed for beeswarm and boxplot)
clear var condition_PHYS condition_EGTA condition_EGTAK
[condition_PHYS{1:size(interclusterdist_phys_GB,1)}]   = deal('PHYS');
[condition_EGTA{1:size(interclusterdist_egta_GB,1)}]   = deal('EGTA');
[condition_EGTAK{1:size(interclusterdist_egtak_GB,1)}] = deal('EGTA/K+');
conditions = [condition_PHYS condition_EGTA condition_EGTAK];

% Plot intercluster distances estimated from Ripley's curves
fig13 = figure;
subplot(121) % beeswarm plot
h = plotSpread({interclusterdist_phys_GB,interclusterdist_egta_GB,interclusterdist_egtak_GB},...
    'xNames',{'PHYS','EGTA','EGTA/K+'});
%set(h{1},'color','k','markersize',10);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','a-synuclein and VAMP2'});
set(gca,'fontsize',14);

subplot(122) % boxplot
boxplot([interclusterdist_phys_GB; interclusterdist_egta_GB; interclusterdist_egtak_GB],conditions);
ylabel('Intercluster distance (nm)');
title({'Intercluster distance','a-synuclein and VAMP2'});
set(gca,'fontsize',14);

if flagprint
    savefig(fig13,fullfile(path_figures_ripley,'intercluster_distance_a-synuclein_VAMP2.fig'))
end

% One-way ANOVA
p = anova1([interclusterdist_phys_GB; interclusterdist_egta_GB; interclusterdist_egtak_GB],conditions,'off');
if p < 0.05
    disp(['Significant difference in a-synuclein/VAMP2 intercluster distance between conditions (p-value = ' num2str(p) ')'])
else
    disp(['No significant difference in a-synuclein/VAMP2 intercluster distance between conditions (p-value = ' num2str(p) ')'])
end

%% Plot the Ripley's functions

path_ripley_curves = fullfile(path_figures_ripley,'curves');
mkdir(path_ripley_curves);

ymin = 0;
ymax = 6000;

% mCLING


fig1 = RipleyPlot(r_hist,H_all_phys_RC(:,all(~isnan(H_all_phys_RC))),'red','mCling PHYS',[ymin ymax]);
fig2 = RipleyPlot(r_hist,H_all_egta_RC(:,all(~isnan(H_all_egta_RC))),'red','mCling EGTA',[ymin ymax]);
fig3 = RipleyPlot(r_hist,H_all_egtak_RC(:,all(~isnan(H_all_egtak_RC))),'red','mCling EGTA/K+',[ymin ymax]);
fig4 = RipleyPlotMeansConditions(r_hist, H_all_phys_RC(:,all(~isnan(H_all_phys_RC))),...
                                         H_all_egta_RC(:,all(~isnan(H_all_egta_RC))),...
                                         H_all_egtak_RC(:,all(~isnan(H_all_egtak_RC))),...
                                         {'PHYS','EGTA','EGTA/K+'},'Mean of Ripley''s - mCling');

% a-synuclein
fig5 = RipleyPlot(r_hist,H_all_phys_GC(:,all(~isnan(H_all_phys_GC))), 'green','a-synuclein PHYS',[ymin ymax]);
fig6 = RipleyPlot(r_hist,H_all_egta_GC(:,all(~isnan(H_all_egta_GC))), 'green','a-synuclein EGTA',[ymin ymax]);
fig7 = RipleyPlot(r_hist,H_all_egtak_GC(:,all(~isnan(H_all_egtak_GC))),'green','a-synuclein EGTA/K+',[ymin ymax]);
fig8 = RipleyPlotMeansConditions(r_hist, H_all_phys_GC(:,all(~isnan(H_all_phys_GC))),...
                                         H_all_egta_GC(:,all(~isnan(H_all_egta_GC))),...
                                         H_all_egtak_GC(:,all(~isnan(H_all_egtak_GC))),...
                                         {'PHYS','EGTA','EGTA/K+'},'Mean of Ripley''s - a-synuclein');

% VAMP2
fig9  = RipleyPlot(r_hist,H_all_phys_BC(:,all(~isnan(H_all_phys_BC))), 'blue','VAMP2 PHYS',[ymin ymax]);
fig10 = RipleyPlot(r_hist,H_all_egta_BC(:,all(~isnan(H_all_egta_BC))), 'blue','VAMP2 EGTA',[ymin ymax]);
fig11 = RipleyPlot(r_hist,H_all_egtak_BC(:,all(~isnan(H_all_egtak_BC))),'blue','VAMP2 EGTA/K+',[ymin ymax]);
fig12 = RipleyPlotMeansConditions(r_hist, H_all_phys_BC(:,all(~isnan(H_all_phys_BC))),...
                                          H_all_egta_BC(:,all(~isnan(H_all_egta_BC))),...
                                          H_all_egtak_BC(:,all(~isnan(H_all_egtak_BC))),...
                                          {'PHYS','EGTA','EGTA/K+'},'Mean of Ripley''s - VAMP2');

if flagprint
    savefig(fig1 ,fullfile(path_ripley_curves,'ripley_mcling_phys.fig'))
    savefig(fig2 ,fullfile(path_ripley_curves,'ripley_mcling_egta.fig'))
    savefig(fig3 ,fullfile(path_ripley_curves,'ripley_mcling_egtak.fig'))
    savefig(fig4 ,fullfile(path_ripley_curves,'ripley_mcling.fig'))
    savefig(fig5 ,fullfile(path_ripley_curves,'ripley_a-synuclein_phys.fig'))
    savefig(fig6 ,fullfile(path_ripley_curves,'ripley_a-synuclein_egta.fig'))
    savefig(fig7 ,fullfile(path_ripley_curves,'ripley_a-synuclein_egtak.fig'))
    savefig(fig8 ,fullfile(path_ripley_curves,'ripley_a-synuclein.fig'))
    savefig(fig9 ,fullfile(path_ripley_curves,'ripley_VAMP2_phys.fig'))
    savefig(fig10,fullfile(path_ripley_curves,'ripley_VAMP2_egta.fig'))
    savefig(fig11,fullfile(path_ripley_curves,'ripley_VAMP2_egtak.fig'))
    savefig(fig12,fullfile(path_ripley_curves,'ripley_VAMP2.fig'))
end

diary off