function plotsweepresults(varargin)
% PLOTSWEEPRESULTS Plot the results of a sweep performed using "sweepfdfits".
%
% SYNTAX:
% plotsweepresults(sweepResults)
%
% plotsweepresults(figureHandle, ...)
%
% Creates a figure with all the fit parameter values as a function of boundary
% position.
%
% NOTE:
% Can only be used for a one-boundary sweep (1D sweep).
%
% INPUT:
% sweepResults = result of a call to the "sweepfdfits" function.
% figureHandle = optional handle to a figure window. If given, this window
%       will be used to create the plot.
%
% SEE ALSO:
% sweepfdfits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

if isempty(varargin)
    error('Missing argument: no sweepResults given');
end

if ishghandle(varargin{1})
    figure(varargin{1});
    varargin(1) = [];
else
    figure();
end

sweepResults = varargin{1};

switch sweepResults.options.sweepType
    case {'left','right'}
        % ok
    otherwise
        error('plotsweepresults can only be used on the results of a one-boundary sweep.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make plots

paramNames = fieldnames(sweepResults.params);
nPlots = length(paramNames)-1;
nRows = floor(sqrt(nPlots));
nCols = ceil(nPlots / nRows);

iPlot = 1;
for iParam = 1:length(paramNames)
    if ~strcmpi(paramNames{iParam}, 'all')
        subplot(nRows, nCols, iPlot);
        makeParamPlot(sweepResults.params.(paramNames{iParam}), paramNames{iParam});
        iPlot = iPlot + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function makeParamPlot(vals, name)
        if strcmp(sweepResults.options.sweepType, 'left')
            sweepAxis = sweepResults.options.sweepPosLeft;
            xLabelType = 'left';
        else
            sweepAxis = sweepResults.options.sweepPosRight;
            xLabelType = 'right';
        end
        plot(sweepAxis, vals, '.b');
        switch sweepResults.options.sweepBoundaryType
            case 'd'
                xlabel(['d_{' xLabelType '} ({\mu}m)']);
            case 'f'
                xlabel(['F_{' xLabelType '} (pN)']);
        end
        title(name);
    end

end
