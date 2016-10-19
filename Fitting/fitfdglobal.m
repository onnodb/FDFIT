function [fitRes, diagnosticInfo] = fitfdglobal(fdc, varargin)
% FITFDGLOBAL Perform a global fit of the "odijk_d0_f0" model on a set of F,d data.
%
% This function uses 'lsqcurvefit' to perform a global fit to a set of DNA
% force-extension curves. This means that all curves are fitted with the
% Odijk eWLC model *simultaneously*, while sharing the physical parameters
% "Lp" (persistence length) and "S" (stretch modulus).
%
% The model that is actually used for fitting, is the "odijk-d0-f0" model
% (see "fitfd"). This is the Odijk eWLC, inversed (expressing force as a
% function of distance, c.f. the paper), with a distance and force offset
% included per-curve.
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
% Lc = value of Lc (contour length) to use while fitting (should be the theoretical,
%       expected value of Lc) (default: 16.5 um).
% startParams = vector [Lp S d0 F0] with starting values for the fit.
%       (Note: specifying d0 and F0 is optional).
% lBounds = lower bounds for the fit parameters; vector like 'startParams'.
% uBounds = upper bounds for the fit parameters; vector like 'startParams'.
%
% FLAG ARGUMENTS:
% computeConfInt = if given, also computes 95% confidence intervals for the
%       fit parameters. WARNING: can be (extremely) slow for large numbers
%       (>= 20) of datasets.
% noTrim = normally, data points above 30 pN ('trimForce' constant in the
%       code below) are automatically discarded; pass in this flag to
%       disable this behavior.
%
% OUTPUT:
% fitRes = structure with the fit parameters found, as well as the value of
%       'Lc' (contour length) used.
%       If the flag argument 'computeConfInt' was given, then each value has
%       the form [fit lbound ubound], where "lbound" and "ubound" are the
%       bounds of the 95% confidence intervals for each fit parameter.
%   .Lp = persistence length (in nm).
%   .Lc = contour length (in um) (not a fit parameter; was fixed to the value
%           'Lc' as passed into this function; see "Key-Value Pair Arguments").
%   .S = stretching modulus (in pN).
%   .d0 = distance offset (in um).
%   .F0 = force offset (in pN).
% diagnosticInfo = structure with extra diagnostic information about the fit,
%       as returned by 'lsqcurvefit'.
%
% NOTE:
% Fitting other models than "odijk_d0_f0" is not currently possible. Given that
% this model gives the best fits to F,d curves, this shouldn't be an issue.
%
% SEE ALSO:
% lsqcurvefit, fitfd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

defaultArgs = struct(...
                      'computeConfInt',             false ...
                    , 'lBounds',                    [] ...
                    , 'Lc',                         16.5 ...
                    , 'startParams',                [50 1500 0 0] ... % Lp S d0 F0
                    , 'uBounds',                    [] ...
                    , 'noTrim',                     false ...
                    );

args = parseArgs(varargin, defaultArgs, {'computeConfInt','noTrim'});

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

args.startParams = checkParamArg(args.startParams, 0, 'startParams');
args.lBounds     = checkParamArg(args.lBounds, -Inf, 'lBounds');
args.uBounds     = checkParamArg(args.uBounds, +Inf, 'uBounds');

if isempty(args.startParams)
    error('fitfdglobal:MissingArgument', 'Missing argument "startParams": no starting values given for fit parameters.');
end

trimForce = 30;     % Odijk model is valid up to 30 pN; trim forces above that


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do global fit
% The way we perform the 'global fit' works as follows:
%  - We concatenate all data points into one big dataset: the variables
%    'fit_data_d' and 'fit_data_f'.
%  - Additionally, we keep track of the data set that each of these data points
%    originally comes from ('fit_data_i').
%  - We now fit a function to this concatenated data, with (2+2N) parameters.
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

if ~args.noTrim
    disp('Only using data up to 30 pN (use the "noTrim" flag to avoid this).');
end

for i = 1:nDatasets
    cur_data_d = fdc.items{i}.d(:)';
    cur_data_f = fdc.items{i}.f(:)';

    if ~args.noTrim
        cur_data_d = cur_data_d(cur_data_f <= trimForce);
        cur_data_f = cur_data_f(cur_data_f <= trimForce);
    end

    fit_data_d = [fit_data_d cur_data_d]; %#ok
    fit_data_f = [fit_data_f cur_data_f]; %#ok
    fit_data_i = [fit_data_i i .* ones(1,length(cur_data_d))]; %#ok
end

nTotalPoints = length(fit_data_d);
assert(length(fit_data_f) == nTotalPoints);
assert(length(fit_data_i) == nTotalPoints);


% Set up the fit.

p0 = [args.startParams(1:2) repmat(args.startParams(3:4), 1, nDatasets)];
%    [Lp S d0(1) F0(1) d0(2) F0(2) ... d0(N) F0 (N)]
%              ^ dataset index

if isempty(args.lBounds)
    lb = [];
else
    lb = [args.lBounds(1:2) repmat(args.lBounds(3:4), 1, nDatasets)];
end
if isempty(args.uBounds)
    ub = [];
else
    ub = [args.uBounds(1:2) repmat(args.uBounds(3:4), 1, nDatasets)];
end

% Set up Jacobian pattern
jp_points = [];
cumTotalNDataPoints = 0;
for i = 1:nDatasets
    nPointsInCurDataSet = sum(fit_data_i == i);
    for j = 1:nPointsInCurDataSet
        jp_points = [jp_points                      ; ...
                     cumTotalNDataPoints+j, 1       ; ...
                     cumTotalNDataPoints+j, 2       ; ...
                     cumTotalNDataPoints+j, 2*i+1   ; ...
                     cumTotalNDataPoints+j, 2*i+2     ...
                     ]; %#ok
    end
    cumTotalNDataPoints = cumTotalNDataPoints + nPointsInCurDataSet;
end
jacobpattern = sparse(jp_points(:,1), jp_points(:,2), ones(length(jp_points),1));

fitopts = optimoptions('lsqcurvefit' ...
                , 'Algorithm',              'trust-region-reflective' ...
                , 'JacobPattern',           jacobpattern ...
                , 'MaxFunEvals',            +Inf ...
                );

[fitParams, ~, residual, exitFlag, output, ~, jacobian] = ...
                lsqcurvefit(...
                    @(x, xdata) fitFun(x, xdata, fit_data_i), ...
                    p0, ...
                    fit_data_d, ...
                    fit_data_f, ...
                    lb, ub, ...
                    fitopts ...
                    );

if args.computeConfInt
    disp('Computing confidence intervals...');
    ci = nlparci(fitParams, residual, 'jacobian', jacobian);

    confInt = struct(...
                     'Lp',              ci(1,:), ...
                     'S',               ci(2,:), ...
                     'd0',              ci(3:2:end-1,:), ...
                     'F0',              ci(4:2:end,:) ...
                     );
else
    confInt = struct(...
                     'Lp',              [], ...
                     'S',               [], ...
                     'd0',              [], ...
                     'F0',              [] ...
                     );
end

fitRes = struct(...
                'Lp',               [fitParams(1) confInt.Lp], ...
                'Lc',               args.Lc, ...
                'S',                [fitParams(2) confInt.S], ...
                'd0',               [fitParams(3:2:end-1)' confInt.d0], ...
                'F0',               [fitParams(4:2:end)' confInt.F0] ...
                );

diagnosticInfo = struct('exitFlag', exitFlag, 'output', output);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [F] = fitFun(p, d, data_i)
        F = zeros(1,nTotalPoints);
        for k = 1:nDatasets
            F(data_i==k) = fOdijkInv_d0_f0(d(data_i==k), p(1), args.Lc, p(2), p(k*2+1), p(k*2+2));
        end
    end

    function [arg] = checkParamArg(arg, autoFillValue, argName)
        switch length(arg)
            case 2
                arg = [arg autoFillValue autoFillValue];
            case 3
                arg = [arg autoFillValue];
            case {0,4}
                % ok
            otherwise
                error('fitfdglobal:InvalidArgument', 'Invalid argument "%s": vector [Lp S d0 F0] expected.', argName);
        end
    end

end
