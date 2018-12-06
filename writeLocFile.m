% This function writes away a table containing localisation data as a
% localisation file using a specified format.
%
% INPUT:
% 	locs..: a table containing localizations (as read in by the
%                function readLocFile.m)
% 	path..: full path to output file
% 	format: 'thunderstorm'
% 	dim...: (optional) dimension (2 or 3 for 2D or 3D localization
%                 files, 2D by default)
%
% Ideas for improvement:
%    - add other formats
%    - add '-1', '-2', ... to filename if name already exists
%
% Author: Ezra Bruggeman, Laser Analytics Group
%
% Last updated on 2 Sept 2018

function [] = writeLocFile(locs,path,format,dim)

if (nargin<4 || isempty(dim)), dim = 2; end

if strcmp(format,'thunderstorm')
    fid = fopen(path,'wt');
    if dim == 2
        header = 'frame,x [nm],y [nm],sigma [nm],intensity [photon],offset [photon],bkgstd [photon],uncertainty [nm]\n';
        fprintf(fid,header);
        if fid > 0
            fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f\n',table2array(locs)');
            fclose(fid);
        end
    elseif dim == 3
        header = 'frame,x [nm],y [nm],z [nm],sigma [nm],intensity [photon],offset [photon],bkgstd [photon],uncertainty [nm]\n';
        fprintf(fid,header);
        if fid > 0
            fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f,%f\n',table2array(locs)');
            fclose(fid);
        end
    end
    
elseif strcmp(format,'rapidstorm')
    disp('Writing away RapidSTORM files is not yet implemented, sorry!')
    
else
    disp('Please select a valid format in writeLocFile.')
    return
end
end