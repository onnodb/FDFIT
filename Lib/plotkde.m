function plotkde(varargin)
% Plots a Kernel Density Estimate of a dataset.
%
% Calls the function 'ksdensity' (from the Statistics Toolbox) to generate
% a kernel smoothing estimate of the data, and generates a nice plot that
% includes the data points themselves, as well.
%
% SYNTAX:
% plotkde(y);
%   Plot the Kernel Density Estimate of a dataset 'y' (a numeric vector).
%
% plotkde(..., 'ksdensityParams', {...});
%   Pass a set of parameters into the 'ksdensity' function. The parameters
%   are given as a cell array.
%
% plotkde(axis_handle, ...);
%   Plot the graph in a specific set of axes.
%
% INPUT:
% y = vector with data points
%
% KEY-VALUE PAIR ARGUMENTS:
% color = color of the plot, as a vector [R G B]. (Default: [0 0 1], blue)
% ksdensityParams = optional set of parameters that should be passed into the
%       'ksdensity' function. Cell array.
% pointPlotScalingFactor = the individual points are drawn near the bottom
%       of the plot, at a height of 'MaxYlim/pointPlotScalingFactor', with
%       'MaxYLim' being the height of the plot ordinate axis.
%       (Default: 20)
%
% FLAG ARGUMENTS:
% plotPoints = if given (default), plots the individual data points as well.
%
% SEE ALSO:
% ksdensity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process & validate input

if isempty(varargin)
    error('Missing arguments: No arguments given');
end

% Explicit axes handle given?
if ishandle(varargin{1})
    ax = varargin{1};
    varargin = varargin(2:end);
else
    ax = gca();
end

% Pop off the first (remaining) argument, that should be the data.
if isempty(varargin)
    error('Missing argument: Data missing');
end
y = varargin{1};
varargin = varargin(2:end);

% Parse remaining key-value pair arguments
defArgs = struct(...
                  'color',                              [0 0 1] ...
                , 'fitNormal',                          false ...
                , 'ksdensityParams',                    {{}} ...
                , 'plotPoints',                         true ...
                , 'pointPlotScalingFactor',             20 ...
                );
args = parseArgs(varargin, defArgs, {'fitNormal', 'plotPoints'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make plot

[f,xi] = ksdensity(y, args.ksdensityParams{:});

plot(ax, xi, f, '-', 'color', args.color);

if args.fitNormal
    pd = fitdist(y, 'Normal');

    hold('on');
    plot(ax, xi, pdf(pd, xi), 'color', [1 0 0]);

    fprintf('Fitted Normal distribution:\n');
    fprintf('  mu    = %g\n', pd.mu);
    fprintf('  sigma = %g\n', pd.sigma);
end

ylims = ylim(ax);

if args.plotPoints
    hold(ax, 'on');
    plot(ax, y, ylims(2)/args.pointPlotScalingFactor*ones(size(y)), '.', 'color', args.color);
end

end
