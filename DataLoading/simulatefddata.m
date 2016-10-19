function [fd] = simulatefddata(varargin)
% SIMULATEFDDATA Generates simulated F,d data by using a model function, and adding some Gaussian noise
%
% SYNTAX:
% fd = simulatefddata(..., 'key', 'value, ...)
%
% KEY-VALUE PAIR ARGUMENTS:
% model = 'odijk'|'twlc'|function handle
%       If a function handle, the function should have the prototype:
%           y = fun(x, param1, param2...)
%       How (x,y) map to (d,F) is set by 'modelType'.
% modelType = 'd(F)'|'F(d)' (only needed if 'model' is a function handle)
% modelParams = (cell) array with parameters for the model function.
% fRange = vector [min max], defining the range of force values.
% dRange = vector [min max], defining the range of distance values.
% nPoints = number of points to generate; note that the *actual* number of
%       points in the simulated curve may be lower (due to removal of invalid
%       and out-of-range points). (Default: 500)
% noise = vector [F_noise, d_noise], with the Gaussian noise amplitude to be
%       added to the force (in pN) and distance values (in um), respectively.
%       'd_noise' is optional, and defaults to 0.
%       Default: [0.2 0.005] (realistic values for optical tweezers
%       experiments).
% offset = vector [F_offset, d_offset]. If non-empty and non-zero, a random
%       offset is added to F and d, respectively.The offsets are taken from
%       a Gaussian distribution with the given amplitudes. (Default: empty).
%       (A good value for 'typical data' could be [0.2 0.05]).
%
% EXAMPLES:
% fd = simulatefddata();
%   This generates a model curve based on the Odijk model, with
%   (Lp;Lc;S) = (50 nm; 16.5 um; 1500 pN), up to 30 pN force.
% fd = simulatefddata('model', 'twlc')
%   Generates a model curve based on the tWLC model, up to 65 pN force.
% fd = simulatefddata('model', @myModelFun, 'modelType', 'F(d)', 'modelParams', [...])
%   Uses a custom model function to generate the data.
%
% NOTES:
% For models of type 'd(F)', interpolation is done in order to get a nice,
% evenly spaced distance axis. This interpolation may add a small (additional)
% amount of noise to the generated data.
%
% SEE ALSO:
% readtwomdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

defArgs = struct(...
                  'dRange',                             [] ...
                , 'fRange',                             [] ...
                , 'model',                              'odijk' ...
                , 'modelParams',                        [] ...
                , 'modelType',                          [] ...
                , 'nPoints',                            500 ...
                , 'noise',                              [0.2 0.005] ...
                , 'offset',                             [] ...
                );
args = parseArgs(varargin, defArgs);

if ischar(args.model)
    switch args.model
        case 'odijk'
            args.modelType = 'd(F)';
            if isempty(args.modelParams)
                args.modelParams = {50 16.5 1500};
            end
            if isempty(args.dRange)
                args.dRange = [10 17.5];
            end
            if isempty(args.fRange)
                args.fRange = [0 30];
            end
            args.model = @fOdijk;
        case 'twlc'
            args.modelType = 'd(F)';
            if isempty(args.modelParams)
                args.modelParams = {50 16.5 1500 440 -637 17 30.6};
            elseif numel(args.modelParams) == 3
                if iscell(args.modelParams)
                    args.modelParams = [args.modelParams {440 -637 17 30.6}];
                else
                    args.modelParams = [args.modelParams 440 -637 17 30.6];
                end
            end
            if isempty(args.dRange)
                args.dRange = [10 17.5];
            end
            if isempty(args.fRange)
                args.fRange = [0 65];
            end
            args.model = @ftWLC;
        otherwise
            error('Invalid argument "model": %s', args.model);
    end
elseif isa(args.model, 'function_handle')
    switch lower(args.modelType)
        case {'d(f)','f(d)'}
            % ok
        otherwise
            if isempty(args.modelType)
                error('Missing argument "modelType"');
            else
                error('Invalid argument "modelType": %s', args.modelType);
            end
    end
else
    error('Invalid argument "model": string or function handle expected');
end

if isempty(args.modelParams)
    error('Missing argument "modelParams"');
elseif isnumeric(args.modelParams)
    % Auto-convert to a cell array
    args.modelParams = num2cell(args.modelParams);
end

if ~isempty(args.dRange) && numel(args.dRange) ~= 2
    error('Invalid argument "dRange": 2-element vector [min max] expected');
end

if ~isempty(args.fRange) && numel(args.fRange) ~= 2
    error('Invalid argument "fRange": 2-element vector [min max] expected');
end

if isempty(args.noise)
    args.noise = [0 0];
elseif isnumeric(args.noise) && numel(args.noise) == 1
    args.noise = [args.noise 0];
elseif isnumeric(args.noise) && numel(args.noise) == 2
    % ok
else
    error('Invalid argument "noise": 2-element vector [noiseF noiseD] expected');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data from model function

metaData = struct(...
                  'file',               '<SimulatedData>' ...
                , 'model',              func2str(args.model) ...
                , 'modelType',          args.modelType ...
                , 'modelParams',        cell2mat(args.modelParams) ...
                , 'dRange',             args.dRange ...
                , 'fRange',             args.fRange ...
                , 'nPoints',            args.nPoints ...
                , 'noise',              args.noise ...
                , 'offset',             args.offset ...
                );

data_d = linspace(args.dRange(1), args.dRange(2), args.nPoints);

switch lower(args.modelType)
    case 'f(d)'
        data_f = args.model(data_d, args.modelParams{:});

    case 'd(f)'
        tmp_f = linspace(args.fRange(1), args.fRange(2), args.nPoints);
        tmp_d = args.model(tmp_f, args.modelParams{:});

        % Remove invalid points before interpolation
        tmp_f(~(isreal(tmp_d) & isfinite(tmp_d))) = [];
        tmp_d(~(isreal(tmp_d) & isfinite(tmp_d))) = [];

        % Use spline interpolation to get evenly-spaced distance axis
        p = spline(tmp_d, tmp_f);
        data_f = ppval(p, data_d);

    otherwise
        error('Invalid argument "modelType": %s', args.modelType);
end

% Remove any out-of-range/NaN/complex values
pointsToRemove = (data_d < args.dRange(1)) | (data_d > args.dRange(2)) | ...
                 (data_f < args.fRange(1)) | (data_f > args.fRange(2)) | ...
                 isnan(data_d) | (~isreal(data_d)) | ...
                 isnan(data_f) | (~isreal(data_f)) ...
                 ;
data_d(pointsToRemove) = [];
data_f(pointsToRemove) = [];

% Add (dummy) time axis
data_t = 1:length(data_d);

% Add noise
data_f = data_f + randn(size(data_f))*args.noise(1);
data_d = data_d + randn(size(data_d))*args.noise(2);

% Add offset
if ~isempty(args.offset)
    data_f = data_f + randn()*args.offset(1);
    data_d = data_d + randn()*args.offset(2);
end

% Make final object
fd = FdData('Simulated data', data_f, data_d, data_t, metaData);


end
