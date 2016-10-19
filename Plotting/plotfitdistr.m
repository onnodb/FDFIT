function plotfitdistr(fits, varargin)
% PLOTFITDISTR Plot distributions of fit values.
%
% SYNTAX:
% plotfitdistr(fits);
% plotfitdistr(..., 'key', 'value', ...);
%
% INPUT:
% fits = results from a call to 'fitfd' on an FdDataCollection.
%
% KEY-VALUE PAIR ARGUMENTS:
% type = hist|kde
%       If 'hist', plots histograms.
%       If 'kde', uses the 'plotkde' function to plot Kernel Density Estimates.
% nBins = number of bins in the histogram (default: 10)
% restrictData = if given, can be used to restrict data to a certain range.
%       (Useful for removing obvious outliers).
%       Should be a cell array, with the form:
%           { {'Lp', [10 50]}, {'paramName', [min max]}, ... }
%
% FLAG ARGUMENTS:
% fitNormal = fits normal distributions to the histograms/kde plots.
% useSubPlots = combine plots into one figure window with subplots.
%
% SEE ALSO:
% fitfd, plotkde

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

defArgs = struct(...
                  'fitNormal',                      false ...
                , 'nBins',                          10 ...
                , 'restrictData',                   [] ...
                , 'type',                           'hist' ...
                , 'useSubPlots',                    false ...
                );
args = parseArgs(varargin, defArgs, {'fitNormal', 'useSubPlots'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect parameter names; we plot one graph per parameter

if isfield(fits, 'fitobject')
    paramNames = coeffnames(fits(1).fitobject);
else
    % Result of a call to 'fitfd' with the 'skipFullFitInfo' flag.
    paramNames = fieldnames(fits);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make plots

if args.useSubPlots
    figure;
    subPlotRows = floor(sqrt(length(paramNames)));
    subPlotCols = ceil(length(paramNames) / subPlotRows);
end
for iParam = 1:length(paramNames)
    curParam = paramNames{iParam};
    if args.useSubPlots
        subplot(subPlotRows, subPlotCols, iParam);
    else
        figure;
    end
    fprintf('### %s\n', curParam);
    makePlot(curParam);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function makePlot(paramName)
        paramVals = cell2mat({fits.(paramName)});

        % Restrict data, if requested.
        restrictionRange = getRestrictionRange(paramName);
        if ~isempty(restrictionRange)
            lengthBeforeRestriction = length(paramVals);
            paramVals(paramVals < restrictionRange(1)) = [];
            paramVals(paramVals > restrictionRange(2)) = [];
            fprintf('Data for %s restricted: %d point(s) removed\n', ...
                    paramName, lengthBeforeRestriction-length(paramVals));
        end

        switch args.type
            case 'hist'
                % Plot the histogram.
                hist(paramVals, args.nBins);

                % Fit normal distribution.
                if args.fitNormal
                    distr = fitdist(paramVals(:), 'normal');

                    % Print properties of Gaussian distribution to screen.
                    fprintf('Fitted Normal distribution:\n');
                    fprintf('  mu    = %g\n', distr.mu);
                    fprintf('  sigma = %g\n', distr.sigma);

                    % Prepare fit plot.
                    xVals = linspace(min(paramVals), max(paramVals), 100);
                    yVals = pdf(distr, xVals);

                    % Rescale the fit to match the histogram.
                    yScale = ylim;
                    yVals = yVals / max(yVals) * yScale(2);

                    % Plot the Gaussian fit.
                    hold('on');
                    plot(xVals, yVals, '-r');
                    hold('off');
                end
            case 'kde'
                if args.fitNormal
                    plotkdeArgs = {'fitNormal'};
                else
                    plotkdeArgs = {};
                end

                plotkde(paramVals(:), plotkdeArgs{:});
            otherwise
                error('Invalid argument "type"');
        end

        % Labels.
        xlabel(paramName);
        ylabel('count');
        title(paramName);
    end

    function [range] = getRestrictionRange(paramName)
        for i = 1:length(args.restrictData)
            if strcmpi(args.restrictData{i}{1}, paramName)
                range = args.restrictData{i}{2};
                return
            end
        end
        range = {};
    end

end
