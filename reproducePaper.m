function reproducePaper(pathToData, whatToReproduce)
% REPRODUCEPAPER Call this function to reproduce all figures from the paper.
%
% SYNTAX:
% reproducePaper /path/to/hdf5/data main
%   Reproduce only the main figures.
% reproducePaper /path/to/hdf5/data sm
%   Reproduce only the Supplementary figures.
% reproducePaper /path/to/hdf5/data all
%   Reproduce all figures.
%
% INPUT:
% pathToData = folder containing all HDF5 files from the data publication
%       referenced in the paper.
% whatToReproduce = main|si|all (default: all)

if nargin < 1 || isempty(pathToData) || ~exist(pathToData, 'dir')
    fprintf('Please select one of the H5 files from the data publication\n');
    fprintf('in the dialog that follows. Alternatively, you can call\n');
    fprintf('this function like so:\n\n');
    fprintf('   >> reproducePaper %s\n', fullfile('path','to','H5','data'));
    fprintf('\nPress [Enter] to continue... ');
    input('', 's');

    [~, pathName, ~] = uigetfile({'*.h5', 'HDF5 files'}, 'Select first data file');
    [pathToData, ~, ~] = fileparts(pathName);
end

if nargin < 2 || isempty(whatToReproduce)
    whatToReproduce = 'all';
else
    whatToReproduce = lower(whatToReproduce);
    switch whatToReproduce
        case {'main','sm','all'}
            % ok
        otherwise
            warning('Unknown item "%s"; assuming "all".', whatToReproduce);
            whatToReproduce = 'all';
    end
end

% Read data.
if ~exist(pathToData, 'dir')
    error(['Path not found. Please specify the path to the folder containing ' ...
           'the H5 files from the data publication referenced in the paper.']);
end

fprintf('Loading data...\n');
fdc = FdDataCollection();
files = dir(fullfile(pathToData, '*.h5'));
if isempty(files)
    error(['No H5 files found. Please specify the path to the folder containing ' ...
           'the H5 files from the data publication referenced in the paper.']);
end
for i = 1:length(files)
    fprintf('.');
    try
        fd = importfddata(fullfile(pathToData, files(i).name));
    catch err
        fprintf('\n!!! Error loading data file %s\n', files(i).name);
        disp(getReport(err, 'full'));
        fprintf('Please check if you have selected the correct data folder, \n');
        fprintf('and if the files have been downloaded correctly.\n');
        error('Error loading data file %s', files(i).name);
    end

    fdc.add(fd);
end
fprintf('\n');


% Make figures.
makeFigures(fdc, ['-' whatToReproduce]);

end
