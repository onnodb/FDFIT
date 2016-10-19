function [mgConc] = fdTagToMgConc(fd)
% FDTAGTOMGCONC Find the associated magnesium concentration for an FdData object.
%
% The force-extension data used in the paper is contained in a series of
% FdData objects that all have their magnesium concentration encoded in
% the object's "tags" property. This function finds such a magnesium-concentration-
% encoding tag in an FdData object's "tags" property, and returns the
% corresponding magnesium concentration (in mM).
%
% SYNTAX:
% mgConc = fdTagToMgConc(fd);
%
% INPUT:
% fd = an FdData object.
%
% OUTPUT:
% mgConc = magnesium concentration (in mM; range 0-999).
%
% NOTE:
% An error is raised if no magnesium concentration tag was found.
%
% SEE ALSO:
% mgConcToFdTag, FdData

for iTag = 1:length(fd.tags)
    if strcmpi(fd.tags{iTag}, 'buffer')
        mgConc = 0;
        return
    elseif (length(fd.tags{iTag}) > 2) && (strcmpi(fd.tags{iTag}(1:2), 'mg'))
        mgConc = sscanf(fd.tags{iTag}(3:end), '%f');
        return
    end
end

error('No magnesium concentration tag found for the fd object "%s".', fd.name);

end

