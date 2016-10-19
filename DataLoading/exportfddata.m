function exportfddata(fd, filename)
% EXPORTFDDATA Persist an FdData object to an HDF5 file
%
% SYNTAX:
% exportfddata(fd, filename)
%
% INPUT:
% fd = an FdData object.
% filename = name of the file to write. NOTE: The 'tilde' (~) is NOT expanded
%       into the current user's home folder on UNIX systems.
%
% SEE ALSO:
% FdData

f = SimpleHdf5File(filename, 'w');

f.writeData('FdtData', [fd.f fd.d fd.t]);

f.setAttributes('/', struct('Name', fd.name));
f.setAttributes('/FdtData', fd.metaData);

writeIfNotEmpty(f, '/Tags', fd.tags);
writeIfNotEmpty(f, '/Marks', fd.marks');
writeIfNotEmpty(f, '/History', fd.history);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function writeIfNotEmpty(f, path, data)
        if isempty(data) || (isstruct(data) && isempty(fieldnames(data)))
            % skip
        else
            f.writeData(path, data);
        end
    end

end
