function [figs] = makeFigures(data, varargin)
% Convenience script for generating figures from figure scripts.
%
% SYNTAX:
% makeFigures()
% makeFigures(data)
% makeFigures(___, figure, ___)
%   Specify figures to generate.
% makeFigures(___, '-main')
%   Generate all main figures.
% makeFigures(___, '-sm')
%   Generate all Supplementary figures.
% makeFigures(___, '-all')
%   Generate all main + Supplementary figures.
% makeFigures(___, 'key', 'value)
%
% INPUT:
% data = either an FdDataCollection, or a TwlcVsMgAnalysis object, or a path
%       to a MAT file containing either. If left empty, an attempt is made
%       to load a variable "data" from the global workspace.
% figure = string representing figure that should be generated (the "XX"
%       part in the "FigureXX.m" filenames).
%
% FLAG ARGUMENTS:
% autoExport = (see 'Figure' class)
%
% OUTPUT:
% figs = struct with created figure objects.
%
% EXAMPLE:
% >> makeFigures(data, '1', 3, 4)
% Generate figures 1, 3 and 4, from "Figure1.m", "Figure3.m" and "Figure4.m",
% respectively.
% >> makeFigures(data, '-all')

validFigures = {{'1','3','4','5'}, ...
                {'DataAlignment','Figure4_S','FmaxSweep','GlobalFitPerformance', ...
                 'GvsMagnesium','JustFitEwlc','ManualGFits','MgVsFmax', ...
                 'MgVsOSPlateauForce'} ...
                };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

%% Data object

if nargin < 1 || isempty(data)
    varNames = {'tvma', 'data', 'c'};
    while isempty(data) && ~isempty(varNames)
        try
            data = evalin('base', varNames{1});
            varNameUsed = varNames{1};
        catch
        end
        varNames(1) = []; % pop
    end

    if isempty(data)
        error('No data object found in global workspace.');
    else
        fprintf('Loaded data object from variable "%s" in global workspace.\n', varNameUsed);
    end
end

if ischar(data)
    % Try and load data object from MAT file.
    fprintf('Loading data object from:\n   %s\n', data);
    data = load(data);

    % Pick the first variable stored in the MAT file.
    fn = fieldnames(data);
    if isempty(fn)
        error('Invalid data.');
    end
    data = data.(fn{1});
end

if isa(data, 'TwlcVsMgAnalysis')
    tvma = data;
    data = tvma.data;
elseif isa(data, 'FdDataCollection')
    tvma = [];
else
    error('Invalid data: FdDataCollection or TwlcVsMgAnalysis object expected.');
end


%% Figure names
figures = {};
while ~isempty(varargin) && isFigureName(varargin{1})
    if strcmpi(varargin{1}, '-main')
        figures = [figures validFigures{1}]; %#ok
    elseif strcmpi(varargin{1}, '-sm')
        figures = [figures validFigures{2}]; %#ok
    elseif strcmpi(varargin{1}, '-all')
        figures = [figures validFigures{1} validFigures{2}]; %#ok
    else
        figures{end+1} = varargin{1}; %#ok
    end
    varargin(1) = []; % pop

    if isnumeric(figures{end})
        figures{end} = num2str(figures{end});
    end
end

if isempty(figures)
    error('Nothing to do: no figures given.');
end


%% Any other arguments (key-value pairs/flags)
defArgs = struct(...
                  'autoExport',                         false ...
                );
args = parseArgs(varargin, defArgs, {'autoExport'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-process data, if needed

if ~data.getByTag('!nooverstretchingdata').isempty()
    disp('Removing data without overstretching plateau.');
    data.remove(data.getByTag('!nooverstretchingdata'));

    if ~isempty(tvma)
        % Will need to re-do analysis...
        warning('Data has changed; re-doing TwlcVsMgAnalysis.');
        tvma = [];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run TwlcVsMgAnalysis, if needed

if isempty(tvma)
    disp('##### Running TwlcVsMgAnalysis...');
    tvma = TwlcVsMgAnalysis('data', data);
    tvma.analyze();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate figures

figs = struct();

commonArgs = {'data', data};
if args.autoExport
    commonArgs = [commonArgs {'autoExport'}];
end

for iFig = 1:length(figures)
    curFig = figures{iFig};

    fprintf('##### Generating Figure %s...\n', curFig);

    fn = ['f' curFig];

    switch lower(curFig)
        case '1'
            figs.(fn) = Figure1(commonArgs{:}, 'tvma', tvma);
        case '3'
            figs.(fn) = Figure3(commonArgs{:}, 'dataAligned', tvma.dataAligned);
        case '4'
            figs.(fn) = Figure4(commonArgs{:});
        case '5'
            figs.(fn) = Figure5(commonArgs{:}, 'tvma', tvma);
        case 'dataalignment'
            figs.(fn) = DataAlignment(commonArgs{:});
        case 'figure4_s'
            figs.(fn) = Figure4_S(commonArgs{:});
        case 'fmaxsweep'
            figs.(fn) = FmaxSweep(commonArgs{:}, 'tvma', tvma);
        case 'globalfitperformance'
            figs.(fn) = GlobalFitPerformance(commonArgs{:});
        case 'gvsmagnesium'
            figs.(fn) = GvsMagnesium(commonArgs{:}, 'tvma', tvma);
        case 'justfitewlc'
            figs.(fn) = JustFitEwlc(commonArgs{:}, 'tvma', tvma);
        case 'manualgfits'
            figs.(fn) = ManualGFits(commonArgs{:}, 'tvma', tvma);
        case 'mgvsfmax'
            figs.(fn) = MgVsFmax(commonArgs{:}, 'tvma', tvma);
        case 'mgvsosplateauforce'
            figs.(fn) = MgVsOSPlateauForce(commonArgs{:});
        otherwise
            error('Unknown figure "%s".', curFig);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [b] = isFigureName(s)
        b = (isnumeric(s) && isscalar(s)) ...
            || ...
            any(strcmpi(s, validFigures{1})) ...
            || ...
            any(strcmpi(s, validFigures{2})) ...
            || ...
            any(strcmpi(s, {'-all','-sm','-main'})) ...
            ;
    end

end
