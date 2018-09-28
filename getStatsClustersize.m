% This function estimates the clustersizes from a bunch of H = L(r) - r
% Ripley's functions. It returns the clustersizes, mean and sd.
% 
% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated on 27 Sept 2018

function [clustersize, Mean, SD] = getStatsClustersize(r_hist,H_all)

% Initialize
clustersize = zeros(size(H_all,2),1);

% Loop over the curves in H_all
for i = 1:size(H_all,2)
    if isnan(H_all(:,i))
        % Sometimes, the H calculated from a region will just be a column
        % filled with NaN. If this is the case, put clustersize = -1
        % (clearly not a real value) to still keep the same dimensions for
        % the output table in which all the different measurements
        % (overlaps, clustersizes etc.) are put together (in the
        % compareConditions.m script). The rows with a clustersize of -1
        % can later be deleted from the table. This way the clustersize can
        % still be linked to the correct synaptosomeID.
        clustersize(i) = -1;
    else
        % Get the clustersize from the position of the max of the curve
        [idx_x,~] = getCoordinatesMax(r_hist,H_all(:,i));
        clustersize(i) = idx_x;
    end
end

% Calculate the mean and sd of the calculated clustersizes
Mean = mean(clustersize);
SD = std(clustersize);

end