function [sweepResults] = sweepfdfits(fd, varargin)
% SWEEPFDFITS Perform fit sweeps on DNA force-extension data.
%
% This function performs a series of fits (using the 'fitfd' function),
% each fit for a different subset of the data. For example, this function
% can be used to cut off different amounts of data in the low-force tail
% of a dsDNA force-extension curve, and study the effect on the fit results.
%
% What subset of the data is fitted is determined by one or two moving
% 'boundaries', which can be either along the force axis, or along the
% distance axis.
%
% SYNTAX:
% sweepFdFits(fd)
% sweepFdFits(..., 'key', value, ...)
%
% INPUT:
% fd = FdData object.
%
% KEY-VALUE PAIR ARGUMENTS:
% model = model to fit; note that this setting is automatically passed
%       into 'fitfd', so there's no need to specify this option in
%       the 'fitfdArgs' argument. See the 'fitfd' function.
% nSweeps = number of sweeps; can be a number, or, for sweepType 'both', a
%       1x2 vector [nSweepsLeft nSweepsRight] (default: 100)
% sweepBoundaries = between which distance/force values the boundaries should
%       be swept. Can be a 1x2 vector, or, for sweepType 'both', a
%       2x2 matrix. See below.
% sweepBoundaryType = 'd' or 'f'; determines whether sweep boundaries are
%       chosen along the distance or the force axis.
% sweepType = 'left'/'right'/'both'; determines which boundary is swept;
%       for 'left'/'right', the other boundary is set to the
%       maximum/minimum of the data. See below. (Default: 'left').
% gofType = which goodness-of-fit statistic to collect. (default: 'rsquare')
% fitfdArgs = cell array with any key-value pair arguments that should be
%       passed into the 'fitfd' function. Optional.
%
% Notes on how the sweep boundaries are determined:
% - For "sweepType" 'both', "sweepBoundaries" is expected to be a 2x2 matrix,
%   with values [startLeft stopLeft; startRight stopRight].
% - For "sweepType" 'left'/'right', "sweepBoundaries" can be either:
%   > a 1x2 matrix specifying the start and stop positions for the moving
%     boundary. The other boundary is set to the leftmost/rightmost data point.
%   > a 2x2 matrix, as above. The appropriate 'start*' position is taken
%     for the static boundary.
%
% OUTPUT:
% A structure with the following fields:
%   .sweeps(iLeft,jRight) = sweep results (iLeft iterates over positions
%                  of the left boundary; etc.)
%       .dataRange
%       .fitExitFlag
%   .params = fit results for the sweeps; each field is a matrix with
%             the same dimensions as .sweeps; values can be NaN
%   .options = the sweep options (properties) actually used, as a
%                    struct; additional fields:
%     .sweepPosLeft = vector of used positions for left boundary
%     .sweepPosRight = vector of used positions for right boundary
%
% SEE ALSO:
% fitfd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

defArgs = struct( ...
                 'model',                       'odijk-f0' ...
               , 'nSweeps',                     100 ...
               , 'sweepBoundaries',             [] ...
               , 'sweepBoundaryType',           'd' ...
               , 'sweepType',                   'left' ...
               , 'gofType',                     'rsquare' ...
               , 'fitfdArgs',                   [] ...
               , 'Lc',                          [] ...
               );

sweepOptions = parseArgs(varargin, defArgs);

if isempty(sweepOptions.fitfdArgs)
    sweepOptions.fitfdArgs = {};
end

switch sweepOptions.sweepType
    case {'left','right'}
        if ~isscalar(sweepOptions.nSweeps)
            error('sweepfdfits:InvalidArgument', 'nSweeps should be a scalar.');
        end
        if ~isequal(size(sweepOptions.sweepBoundaries), [1 2])
            error('sweepfdfits:InvalidArgument', 'sweepBoundaries should be a 1x2 vector.');
        end
    case 'both'
        if ~(isscalar(sweepOptions.nSweeps) || isequal(size(sweepOptions.nSweeps), [1 2]) )
            error('sweepfdfits:InvalidArgument', 'nSweeps should be either a scalar or a 1x2 vector.');
        end
        if ~(isequal(size(sweepOptions.sweepBoundaries), [1 2]) || isequal(size(sweepOptions.sweepBoundaries), [2 2]))
            error('sweepfdfits:InvalidArgument', 'sweepBoundaries should be either a 1x2 vector or a 2x2 matrix.');
        end
    otherwise
        error('sweepfdfits:InvalidArgument', 'Invalid sweepType "%s"; should be "left", "right" or "both".', sweepOptions.sweepType);
end

switch sweepOptions.model
    case 'odijk'
        nParams = 3;
        paramNames = {'Lp','Lc','S'};
    case 'odijk-f0'
        nParams = 4;
        paramNames = {'Lp','Lc','S','F0'};
    case 'odijk-d0'
        nParams = 3;
        paramNames = {'Lp','S','d0'};
    case 'odijk-d0-f0'
        nParams = 4;
        paramNames = {'Lp','S','d0','F0'};
    case 'twlc'
        nParams = 5;
        paramNames = {'Lp','Lc','S','g0','g1'};
    case 'twlc-lc-fixed'
        nParams = 4;
        paramNames = {'Lp','S','g0','g1'};
        if isempty(sweepOptions.Lc)
            error('sweepfdfits:LcMissing', 'Fit requested with Lc fixed, but no Lc value given.');
        end
        sweepOptions.fitfdArgs = [{'Lc', sweepOptions.Lc} sweepOptions.fitfdArgs];
    case 'twlc-f0'
        nParams = 6;
        paramNames = {'Lp','Lc','S','g0','g1','F0'};
end

switch sweepOptions.sweepBoundaryType
    case {'d','f'}
        % ok
    otherwise
        error('sweepfdfits:InvalidArgument', 'Invalid sweepBoundaryType "%s".', sweepOptions.sweepBoundaryType);
end

if isempty(sweepOptions.sweepBoundaries)
    switch sweepOptions.sweepBoundaryType
        case 'd'
            sweepOptions.sweepBoundaries = [min(fd.d) max(fd.d)];
        case 'f'
            sweepOptions.sweepBoundaries = [min(fd.f) max(fd.f)];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up sweeps

% Gather boundary position lists for left and right boundaries; if we're
% not sweeping a boundary, the list will just contain one position.
switch sweepOptions.sweepType
    case 'left'
        sweepPosLeft  = makeSweepPositions('L', sweepOptions.sweepBoundaries, sweepOptions.nSweeps);
        sweepPosRight = makeBoundaryPos('R', sweepOptions.sweepBoundaries);
    case 'right'
        sweepPosLeft  = makeBoundaryPos('L', sweepOptions.sweepBoundaries);
        sweepPosRight = makeSweepPositions('R', sweepOptions.sweepBoundaries, sweepOptions.nSweeps);
    case 'both'
        sweepPosLeft  = makeSweepPositions('L', sweepOptions.sweepBoundaries, sweepOptions.nSweeps);
        sweepPosRight = makeSweepPositions('R', sweepOptions.sweepBoundaries, sweepOptions.nSweeps);
end

% Right boundary should be swept from right to left.
sweepPosRight = fliplr(sweepPosRight);

% Note: the way we produce and gather results is slightly convoluted; this
% is for compatibility with the 'parfor' construct.

% Number of sweeps and eventual size of sweep results matrix.
nSweeps = length(sweepPosLeft)*length(sweepPosRight);
sweepsSize = [length(sweepPosLeft),length(sweepPosRight)];
sweepBoundaryType = sweepOptions.sweepBoundaryType;

% Set up sweep results as 1D vectors / cell arrays.
resSweeps = struct('dataRange', NaN, 'fitExitFlag', NaN);
for iLeft = 1:length(sweepPosLeft)
    for iRight = 1:length(sweepPosRight)
        resSweeps(sub2ind(sweepsSize, iLeft, iRight)).dataRange = [sweepPosLeft(iLeft), sweepPosRight(iRight)];
    end
end

resParams = cell(nSweeps,1);
resErr    = zeros(nSweeps,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Go!

parfor iSweep = 1:nSweeps
    curDataMin = resSweeps(iSweep).dataRange(1);
    curDataMax = resSweeps(iSweep).dataRange(2);

    data_d = [];
    data_f = [];

    switch sweepBoundaryType
        case 'd'
            data_d = fd.d(fd.d >= curDataMin & fd.d <= curDataMax); %#ok
            data_f = fd.f(fd.d >= curDataMin & fd.d <= curDataMax);
        case 'f'
            data_d = fd.d(fd.f >= curDataMin & fd.f <= curDataMax);
            data_f = fd.f(fd.f >= curDataMin & fd.f <= curDataMax);
    end

    f    = [];
    gof  = [];
    fo   = [];

    % Try to fit selected subset of data.
    if length(data_d) > nParams
        % Can only fit if # data points > # free parameters.
        try
            [f, gof, fo] = fitfd(...
                                 FdData('fittmp', data_f, data_d, 1:length(data_f)), ...
                                        'model', sweepOptions.model, ...
                                        sweepOptions.fitfdArgs{:} ...
                                        ); %#ok
        catch err
            if strcmp(err.identifier, 'curvefit:fit:nanComputed') ...
            || strcmp(err.identifier, 'curvefit:fit:complexValueComputed') ...
            || strcmp(err.identifier, 'curvefit:fit:notEnoughDataPoints') ...
            || strcmp(err.identifier, 'fitfd:NotEnoughForceData')
                % not getting a fit; that's OK
            else
                rethrow(err);
            end
        end
    end

    % Collect & save results.
    if isempty(f)
        resSweeps(iSweep).fitExitFlag = NaN;
        resErr(iSweep)                = NaN;
    else
        resSweeps(iSweep).fitExitFlag = fo.exitflag;
        resParams{iSweep}             = coeffvalues(f);
        resErr(iSweep)                = gof.(sweepOptions.gofType);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect results

sweepResults = struct();

sweepResults.sweeps = reshape(resSweeps, sweepsSize);

sweepResults.params.all = reshape(resParams, sweepsSize);
allResParams = cell2mat(resParams);
for iParam = 1:nParams
    sweepResults.params.(paramNames{iParam}) = ...
        reshape(allResParams(:,iParam), sweepsSize);
end
sweepResults.params.err = reshape(resErr, sweepsSize);

% Add options to results struct, too.
sweepOptions.sweepPosLeft  = sweepPosLeft;
sweepOptions.sweepPosRight = sweepPosRight;
sweepResults.options = sweepOptions;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [sweepPos] = makeSweepPositions(boundary, sweepBoundaries, nSweeps)
        if boundary == 'R' && size(sweepBoundaries,1) > 1
            sweepBoundStart = sweepBoundaries(2,1);
            sweepBoundStop  = sweepBoundaries(2,2);
        else
            sweepBoundStart = sweepBoundaries(1,1);
            sweepBoundStop  = sweepBoundaries(1,2);
        end
        if boundary == 'R' && length(nSweeps) > 1
            nSweeps = nSweeps(2);
        else
            nSweeps = nSweeps(1);
        end
        sweepStep = (sweepBoundStop - sweepBoundStart) / (nSweeps-1);
        sweepPos = sweepBoundStart:sweepStep:sweepBoundStop;
        if boundary == 'R'
            sweepPos = sweepPos(end:-1:1);
        end
    end

    function [bPos] = makeBoundaryPos(boundary, sweepBoundaries)
        if boundary == 'R'
            if size(sweepBoundaries,1) > 1
                bPos = sweepBoundaries(2,1);
            else
                switch sweepOptions.sweepBoundaryType
                    case 'd'
                        bPos = max(fd.d);
                    case 'f'
                        bPos = max(fd.f);
                end
            end
        else
            if size(sweepBoundaries,1) > 1
                bPos = sweepBoundaries(1,1);
            else
                switch sweepOptions.sweepBoundaryType
                    case 'd'
                        bPos = min(fd.d);
                    case 'f'
                        bPos = min(fd.f);
                end
            end
        end
    end

end

