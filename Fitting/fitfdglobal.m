function [fitRes, diagnosticInfo] = fitfdglobal(fdc, varargin)
% FITFDGLOBAL Perform a global fit of a fit model on a set of F,d data.
%
% This function uses 'lsqcurvefit' to perform a global fit to a set of DNA
% force-extension curves. This means that all curves are fitted with the
% chosen model *simultaneously*, while sharing the physical parameters
% given by the "sharedParams" argument.
%
% SYNTAX:
% fitRes = fitfdglobal(fdc);
% [fitRes, diagnosticInfo] = fitfdglobal(fdc);
% ... = fitfdglobal(..., 'key', value);
%
% INPUT:
% fdc = an FdDataCollection.
%
% KEY-VALUE PAIR ARGUMENTS:
% sharedParams = cell array of strings, with the names of fit parameters
%       that should be shared among datasets (as opposed to varying per
%       dataset).
%
% FLAG ARGUMENTS:
% computeConfInt = if given, also computes 95% confidence intervals for the
%       fit parameters. WARNING: can be (extremely) slow for large numbers
%       (>= 20) of datasets.
%
% OUTPUT:
% fitRes = structure with the fit parameters found.
%       If the flag argument 'computeConfInt' was given, then each value has
%       the form [fit lbound ubound], where "lbound" and "ubound" are the
%       bounds of the 95% confidence intervals for each fit parameter.
% diagnosticInfo = structure with extra diagnostic information about the fit,
%       as returned by 'lsqcurvefit'.
%
% SEE ALSO:
% lsqcurvefit, fitfd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

defaultArgs = struct(...
                      'computeConfInt',             false ...
                    , 'model',                      [] ...
                    , 'sharedParams',               [] ...
                    );

args = parseArgs(varargin, defaultArgs, {'computeConfInt'});

if isa(fdc, 'FdDataCollection')
    % ok
else
    error('fitfdglobal:InvalidArgument', ...
          'Invalid argument: fdc must be an FdDataCollection.');
end

if fdc.length == 0
    error('fitfdglobal:MissingData', 'No data to fit: FdDataCollection is empty.');
elseif fdc.length == 1
    warning('fitfdglobal:SingleCurve', 'Only fitting a single F,d curve; you may want to use "fitfd" instead.');
end

if isempty(args.model) || ~isa(args.model, 'BasicFdFitModel')
    error('fitfdglobal:InvalidFitModel', 'Invalid fit model.');
end

if isempty(args.sharedParams) || ~iscellstr(args.sharedParams)
    error('fitfdglobal:InvalidSharedParams', 'sharedParams should be a cell array of fit parameter names.');
end
for i = 1:length(args.sharedParams)
    if ~any(strcmp(args.sharedParams{i}, args.model.fitParamNames))
        error('fitfdglobal:InvalidSharedParams', 'sharedParams: invalid fit parameter %s', ...
                args.sharedParams{i});
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do global fit
% The way we perform the 'global fit' works as follows:
%  - We concatenate all data points into one big dataset: the variables
%    'fit_data_d' and 'fit_data_f'.
%  - Additionally, we keep track of the data set that each of these data points
%    originally comes from ('fit_data_i').
%  - We now fit a function to this concatenated data, with (M+N*K) parameters.
%    The first M are the shared parameters
%    The first two of these are "Lp" and "S", which are shared. The remaining
%    2N are "d0" and "F0", which may vary per fitted dataset.
%  - In the fit function "fitFun", we loop over all datasets, and generate
%    a concatenated series of "F" (force) values for the current values of
%    all of the above fit (2+2N) fit parameters.
%  - We pass the structure of the Jacobian into the fit function, which
%    greatly speeds up the fitting process.

nDatasets = fdc.length;

% Concatenate all data.

fit_data_d = [];            % concatenated distance values
fit_data_f = [];            % concatenated force values
fit_data_i = [];            % origin of each data point

for i = 1:nDatasets
    % Check if data is OK for fitting (and possibly filter data).
    fd_fit = args.model.validateData(fdc.items{i});

    cur_data_d = fd_fit.d(:)';
    cur_data_f = fd_fit.f(:)';

    fit_data_d = [fit_data_d cur_data_d]; %#ok
    fit_data_f = [fit_data_f cur_data_f]; %#ok
    fit_data_i = [fit_data_i i .* ones(1,length(cur_data_d))]; %#ok
end

nTotalPoints = length(fit_data_d);
assert(length(fit_data_f) == nTotalPoints);
assert(length(fit_data_i) == nTotalPoints);


% Figure out which parameters to fit, share, and fix.
sharedParams    = args.sharedParams;
nonSharedParams = setdiff(args.model.fitParamNames, args.sharedParams);
nonSharedParams = nonSharedParams(:)';

nSharedParams    = length(sharedParams);
nNonSharedParams = length(nonSharedParams);
nParams          = args.model.nFitParams;
assert(nSharedParams + nNonSharedParams == nParams);

sharedParamMapping = zeros(1,nSharedParams);
for i = 1:nSharedParams
    sharedParamMapping(i) = find(strcmp(sharedParams(i), args.model.fitParamNames));
end

nonSharedParamMapping = zeros(1,nNonSharedParams);
for i = 1:nNonSharedParams
    nonSharedParamMapping(i) = find(strcmp(nonSharedParams(i), args.model.fitParamNames));
end


% Set up the fit.
sharedParamVal = zeros(1,nSharedParams);
sharedParamLow = zeros(1,nSharedParams);
sharedParamUpp = zeros(1,nSharedParams);
for i = 1:nSharedParams
    sharedParamVal(i) = args.model.fitParams.(sharedParams{i});
    sharedParamLow(i) = args.model.fitParamBounds.(sharedParams{i})(1);
    sharedParamUpp(i) = args.model.fitParamBounds.(sharedParams{i})(2);
end

nonSharedParamVal = zeros(1,nNonSharedParams);
nonSharedParamLow = zeros(1,nNonSharedParams);
nonSharedParamUpp = zeros(1,nNonSharedParams);
for i = 1:nNonSharedParams
    nonSharedParamVal(i) = args.model.fitParams.(nonSharedParams{i});
    nonSharedParamLow(i) = args.model.fitParamBounds.(nonSharedParams{i})(1);
    nonSharedParamUpp(i) = args.model.fitParamBounds.(nonSharedParams{i})(2);
end

p0 = [sharedParamVal  repmat(nonSharedParamVal, 1, nDatasets)];
% e.g., for eWLC, with Lp and S shared:
%    [Lp S d0(1) F0(1) d0(2) F0(2) ... d0(N) F0 (N)]
%             ^ dataset index

lb = [sharedParamLow  repmat(nonSharedParamLow, 1, nDatasets)];
ub = [sharedParamUpp  repmat(nonSharedParamUpp, 1, nDatasets)];


% Set up Jacobian pattern
disp('Setting up Jacobian...');
jp_points = zeros(nTotalPoints*nParams, 2);
jp_idx = 1;
for i = 1:nDatasets
    nPointsInCurDataSet = sum(fit_data_i==i);
    for j = 1:nPointsInCurDataSet
        for k = 1:nSharedParams
            jp_points((jp_idx+j-2)*nParams+k,:) = [jp_idx+j-1, k];
        end
        for k = 1:nNonSharedParams
            jp_points((jp_idx+j-2)*nParams+nSharedParams+k,:) = ...
                     [jp_idx+j-1,  nSharedParams+1+(i-1)*nNonSharedParams+(k-1)];
        end
    end
    jp_idx = jp_idx + nPointsInCurDataSet;
end
jacobpattern = sparse(jp_points(:,1), jp_points(:,2), ones(size(jp_points,1),1));

fitopts = optimoptions('lsqcurvefit' ...
                , 'Algorithm',              'trust-region-reflective' ...
                , 'JacobPattern',           jacobpattern ...
                , 'MaxFunEvals',            +Inf ...
                );

% Set up fit function, its parameters, and its data.
switch args.model.dependentVariable
    case 'F'
        depData   = fit_data_f;
        indepData = fit_data_d;
    case 'd'
        depData   = fit_data_d;
        indepData = fit_data_f;
    otherwise
        error('Invalid dependent variable for model object.');
end

fun = args.model.getFitFun();


% Execute fit.
disp('Executing lsqcurvefit...');
[fitParams, ~, residual, exitFlag, output, ~, jacobian] = ...
                lsqcurvefit(...
                    @(x, xdata) fitFun(x, xdata, fun, fit_data_i), ... % f(x)
                    p0, ...                 % start values for fit params
                    indepData, depData, ... % x, y for 'y = f(x)'
                    lb, ub, ...             % lower and upper bounds for fit params
                    fitopts ...
                    );

if args.computeConfInt
    disp('Computing confidence intervals...');
    ci = nlparci(fitParams, residual, 'jacobian', jacobian);

    confInt = paramVectorToStruct(ci);
else
    confInt = [];
end

fitRes = paramVectorToStruct(fitParams);

if ~isempty(confInt)
    fn = fieldnames(fitRes);
    for i = 1:length(fn)
        fitRes.(fn{i}) = [fitRes.(fn{i}), confInt.(fn{i})];
    end
end

diagnosticInfo = struct('exitFlag', exitFlag, 'output', output);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [depVar] = fitFun(paramVals, indepVar, fun, data_i)
        depVar = zeros(size(indepVar));
        for m = 1:nDatasets
            p = cell(1,args.model.nFitParams);
            for q = 1:nSharedParams
                p{sharedParamMapping(q)} = paramVals(q);
            end
            for q = 1:nNonSharedParams
                p{nonSharedParamMapping(q)} = ...
                        paramVals( nSharedParams + nNonSharedParams*(m-1) + q );
            end
            depVar(data_i==m) = fun( p{:}, indepVar(data_i==m) );
        end
    end

    function [res] = paramVectorToStruct(paramVals)
        fn = args.model.fitParamNames;
        res = struct();
        for q = 1:length(fn)
            res.(fn{q}) = [];
        end
        for q = 1:nSharedParams
            res.(fn{sharedParamMapping(q)}) = paramVals( sharedParamMapping(q) );
        end
        for m = 1:nDatasets
            for q = 1:nNonSharedParams
                res.(fn{nonSharedParamMapping(q)}) = ...
                    [ res.(fn{nonSharedParamMapping(q)}), ...
                      paramVals( nSharedParams + nNonSharedParams*(m-1) + q ) ];
            end
        end
    end

end
