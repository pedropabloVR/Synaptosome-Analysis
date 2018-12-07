% This function reads localisation outputfiles from different SMLM
% reconstruction software in as a table.
%
% Files are read in as 2D or 3D depending on the number of columns. For
% ThunderSTORM files, it is assumed that maximum likelihood was used during
% reconstruction (recommended). If least squares was used instead, there
% will be an extra column 'chi2' in the file and it will be read in
% incorrectly (as 3D instead of 2D, or give an error because not enough
% variablenames were specified).
%
% INPUT:
% 	filepath......: full path to localisation file
% 	format........: 'thunderstorm' or 'rapidstorm'
% 	dim (optional): dimension (2 or 3 for 2D or 3D localization files,
%                   2D by default)
%
% OUTPUT:
%   locs..........:
% 
% Ideas for improvement:
%    - add other formats
%
% Author: Ezra Bruggeman, Laser Analytics Group
% 
% Last updated on 8 Sept 2018

function locs = readLocFile(filepath, format, dim)

if (nargin<7 || isempty(dim)), dim = 2; end

if strcmp(format,'thunderstorm')
    if dim == 2
        locs = dlmread(filepath,',',1,0);
        if size(locs,2) == 8
            locs = array2table(locs, 'VariableNames', ...
                {'frame','x','y','sigma','intensity','offset','bkgstd','uncertainty'});
        elseif size(locs,2) == 9
            locs = array2table(locs(:,1:end-1), 'VariableNames', ...
                {'frame','x','y','sigma','intensity','offset','bkgstd','uncertainty'});
        end
    elseif dim == 3
        locs = dlmread(filepath,',',1,0);
        locs = array2table(locs, 'VariableNames',...
            {'frame','x','y','z','sigma','intensity','offset','bkgstd','uncertainty'});
    end
    
elseif strcmp(format,'rapidstorm')
    if dim == 2
        locs = dlmread(filepath,' ',1,0);
        locs = locs(:,1:4);
        locs = array2table(locs, 'VariableNames', ...
            {'x','y','frame','intensity'});
    elseif dim == 3
        locs = dlmread(filepath,' ',1,0);
        locs = locs(:,1:5);
        locs = array2table(locs, 'VariableNames', ...
            {'x','y','z','frame','intensity'});
    end
end
