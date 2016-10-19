function [fd] = importfddata(filename)
% IMPORTFDDATA Import an FdData object from an HDF5 file that was exported using "exportfddata"
%
% SYNTAX:
% fd = importfddata(filename);
%
% SEE ALSO:
% exportfddata

f = SimpleHdf5File(filename, 'r');

fdt = f.readData('FdtData');
metaData = f.getAttributes('/FdtData');
nameData = f.getAttributes('/');
name = nameData.Name;

tags = {};
if f.datasetExists('/Tags')
    tags = f.readData('/Tags');
end

marks = struct();
if f.datasetExists('/Marks')
    marks = f.readData('/Marks');
end

history = {};
if f.datasetExists('/History')
    history = f.readData('/History');
end

fd = FdData(name, fdt(:,1), fdt(:,2), fdt(:,3), metaData, marks, history);

if ~isempty(tags)
    fd.addTag(tags{:});
end

end
