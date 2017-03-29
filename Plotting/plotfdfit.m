function plotfdfit(varargin)
% PLOTFDFIT Plot force-extension data together with an associated fit.
%
% SYNTAX:
% plotfdfit(fd, fitobject, model)
%       Plot the results of a fit made with the "fitfd" function.
% plotfdfit(fd, fitparams, model)
%       Pass a numeric vector with parameter values instead. Useful for fit
%       results as returned from "fitfdglobal".
% plotfdfit(..., 'key', value, ...)
%       For a description of the key-value pair arguments, see below.
% plotfdfit(ax, ...)
%       Plot the results in a specific axes system.
%
% INPUT:
% fd = FdData object
% fitobject = cfit object, as returned by the 'fit' or 'fitfd' functions.
% model = model used for fitting (a BasicFdFitModel descendent).
% params = numeric vector with fit parameters.
% ax = axes handle (optional).
%
% KEY-VALUE PAIR ARGUMENTS:
% style = plot style. Available styles:
%         - 'normal' (default)
%         - 'semilog' (logarithmic F axis)
%         - 'log' (log-log scale)
%
% SEE ALSO:
% fitfd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

if nargin == 0
    error('plotfdfit:InvalidArgument', 'Invalid arguments: no arguments given.');
end

if ishghandle(varargin{1})
    axesHandle = varargin{1};
    varargin(1) = [];
else
    figure();
    axesHandle = gca();
end

if length(varargin) < 3
    error('plotfdfit:InvalidArgument', 'Invalid arguments: no arguments given.');
end

[fd, fitRes, model] = varargin{1:3};
varargin = varargin(4:end);

if isa(fitRes, 'cfit')
    resMode = 'fitobject';
    fitObject = fitRes;
elseif isstruct(fitRes)
    resMode = 'fitparams';
    fitParams = fitRes;
else
    error('plotfdfit:InvalidFitResult', ...
          'Invalid arguments: fit result should be a fit object or a fit struct');
end

if ~isa(model, 'BasicFdFitModel')
    error('plotfdfit:InvalidFitModel', ...
          'Invalid model: object of BasicFdFitModel descendent expected.');
end
switch model.dependentVariable
    case {'F', 'd'}
        % ok
    otherwise
        error('Invalid dependent variable for model object.');
end

% Parse remaining key-value pair arguments
defaultArgs = struct(...
                      'style',              'normal' ...
                    );

args = parseArgs(varargin, defaultArgs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Make plot

switch args.style
    case 'normal'
        plotfun = @plot;
    case 'semilog'
        plotfun = @semilogy;
    case 'log'
        plotfun = @log;
end

% Plot data.
plotfun(axesHandle, fd.d, fd.f, '.b');
hold(axesHandle, 'on');

% Plot fit.
switch resMode

    case 'fitobject'
        switch model.dependentVariable
            case 'F'
                makeFitPlot(fd.d, fitObject(fd.d));
            case 'd'
                makeFitPlot(fitObject(fd.f), fd.f);
        end

        fParamList = [coeffnames(fitObject) ...
                      num2cell(coeffvalues(fitObject))' ...
                      num2cell(confint(fitObject))' ...
                      ]';
        paramText = sprintf('%8.8s: %12g  (%12g -- %-12g)\n', fParamList{:});

    case 'fitparams'
        fun = model.getFitFun();
        p = num2cell(fitResStructToVect(fitRes, model));

        switch(model.dependentVariable)
            case 'F'
                makeFitPlot(fd.d, fun(p{:}, fd.d));
            case 'd'
                makeFitPlot(fun(p{:}, fd.f), fd.f);
        end
        paramText = '';
        fn = fieldnames(fitRes);
        for i = 1:length(fn)
            paramText = [paramText ...
                         sprintf('%8.8s: %12g\n', fn{i}, fitRes.(fn{i}) )]; %#ok
        end
end


text(0.05, 0.85, paramText, ...
        'Parent',      axesHandle, ...
        'Units',       'normalized', ...
        'FontName',    'FixedWidth', ...
        'Interpreter', 'none' ...
        );
hold(axesHandle, 'off');

xlim(axesHandle, [min(fd.d) max(fd.d)]);
ylim(axesHandle, [min(fd.f) max(fd.f)]);
xlabel(axesHandle, 'Distance ({\mu}m)');
ylabel(axesHandle, 'Force (pN)');
if ~isempty(fd.name)
    title(axesHandle, fd.name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function makeFitPlot(d, f)
        % To get a nice line, make sure the data is ordered by ascending 'd'
        sortVals = sortrows([d(:) f(:)]);
        plotfun(axesHandle, sortVals(:,1), sortVals(:,2), '-r');
    end

    function [p] = fitResStructToVect(s, model)
        fn = fieldnames(s);
        p = zeros(1,length(model.fitParamNames));
        for m = 1:length(fn)
            idx = find(strcmp(fn{m}, model.fitParamNames));
            if isempty(idx)
                error('plotfdfit:InvalidArgument', ...
                      'Invalid argument fitparams: unknown model parameter "%s".', fn{m});
            end
            p(idx) = s.(fn{m});
        end
    end

end

