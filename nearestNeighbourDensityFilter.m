% This function performs density-based filtering on a set of points (works
% for 2D or 3D data).
%
% It counts the number of neighbours within a radius 'r' for each point,
% and if it has less than 'n' neighbours, the point is filtered out.
%
% INPUT:
%   X............: nx2 or nx3 double containing coordinates of n points
%   search_radius: radius for which number of neighbours is evaluated
%   N_min........: minimum number of neighbours a point needs to have in
%                  order for it not to be filtered out
%
% OUTPUT:
%   locs_filtered: mx2 or mx3 double containing the coordinates of the
%                  filtered points (with m < or = n, the number of points
%                  from the input)
%   indeces......: nx1 double (with same length as number of points from
%                  the input X) that can be used to check which points were
%                  kept and which were filtered out
%                  (1 = kept, 0 = filtered out)
%   numNeighbours: nx1 double (with same length as number of points from
%                  the input X) that contains the number of neighbours
%                  within a radius search_radius for each localization
%
% Author: Ezra Bruggeman, Laser Analytics Group
%
% Last updated on 5 Sept 2018


function [locs_filtered, indeces, numNeighbours] = ...
    nearestNeighbourDensityFilter(X,search_radius,N_min,show)

disp('Density filtering:')

% Get number of neighbours within a radius r for every localization
[idx,~]= rangesearch(X,X,search_radius);
numNeighbours = cellfun(@(x) numel(x), idx);

% Create an array with the same length as the number of localizations with
% 1 at the index of every localization that has more than n neighbours
indeces = zeros(size(numNeighbours));
indeces(numNeighbours > N_min) = 1;

% Use the mask to filter the localizations
if size(X,2) == 2 % If 2D
    x = X(:,1); y = X(:,2);
    locs = array2table([indeces x y],'VariableNames',{'id','x','y'});
    locs_filtered = locs(locs.id == 1,:);
    locs_filtered = [locs_filtered.x locs_filtered.y];
elseif size(X,2) == 3 % If 3D
    x = X(:,1); y = X(:,2); z = X(:,3);
    locs = array2table([indeces x y z],'VariableNames',{'id','x','y','z'});
    locs_filtered = locs(locs.id == 1,:);
    locs_filtered = [locs_filtered.x locs_filtered.y locs_filtered.z];
end

disp(['    Number of points before filtering: ' num2str(size(X,1))])
disp(['    Number of points after  filtering: ' num2str(size(locs_filtered,1))])
disp(['    Percentage of points filtered out: ' num2str(round(100-(100*size(locs_filtered,1)/size(X,1)))) '%'])


if show
    if size(X,2) == 2
        figure('units','normalized','outerposition',[0 0 1 1]);
        
        xmax = max(locs.x); xmin = min(locs.x);
        ymax = max(locs.y); ymin = min(locs.y);
        
        subplot(131)
        scatter(locs.x,locs.y,'.r')
        title('Before filtering','FontSize',20)
        pbaspect([1 1 1]); xlim([xmin xmax]); ylim([ymin ymax]);
        
        subplot(132)
        scatter(locs_filtered(:,1),locs_filtered(:,2),'.b')
        title('After filtering','FontSize',20)
        pbaspect([1 1 1]); xlim([xmin xmax]); ylim([ymin ymax])
        
        subplot(133)
        scatter(locs.x,locs.y,'.r')
        hold on
        scatter(locs_filtered(:,1),locs_filtered(:,2),'.b')
        title('Overlay','FontSize',20)
        pbaspect([1 1 1]); xlim([xmin xmax]); ylim([ymin ymax])
        
    elseif size(X,2) == 3
        print('No option to display 3D data (yet).')
    end
end

end