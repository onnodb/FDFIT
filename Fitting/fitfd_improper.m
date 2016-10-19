function [varargout] = fitfd_improper(fd, varargin)
% FITFD_IMPROPER Fits DNA force-extension curves using an incorrect approach.
%
% This function is used to fit force-extension data with the Odijk eWLC model,
% but in a representation (distance as a function of force) that results in
% a dependence of the fit results on the range of data used for fitting.
%
% This function should NOT be used for actual data analysis, but is only
% included for comparison. See for example Figure 4.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse input

defaultArgs = struct(...
                      'startParams',        [] ...
                    , 'lBounds',            [] ...
                    , 'uBounds',            [] ...
                    , 'model',              'odijk-d0-f0' ...
                    , 'noTrim',             false ...
                    , 'Lc',                 16.5 ...
                    );
args = parseArgs(varargin, defaultArgs, {'noTrim'});

if ~strcmpi(args.model, 'odijk-d0-f0')
    error('Invalid model "%s": only the "odijk-d0-f0" model is supported by this function.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process & validate input

% Automatically set starting point, if necessary.
if isempty(args.startParams)
    args.startParams = [50 1500 0 0];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do fit

% If fitting Odijk, only use forces up to 30 pN (unless this behavior has been
% disabled explicitly).
if ~args.noTrim
    disp('Only using data up to 30 pN (use the "noTrim" flag to avoid this).');
    fd_fit = fd.subset('f', [-Inf 30]);
else
    fd_fit = fd;
end

% Do actual fit.
ftype = fittype(@(Lp, S, d0, F0, x) fOdijk_d0_f0(x, Lp, args.Lc, S, d0, F0));
fopt = fitoptions(ftype);
fopt.StartPoint = args.startParams;
fopt.Lower      = args.lBounds;
fopt.Upper      = args.uBounds;
fopt.MaxIter    = 120;

[fitobject, gof, output] = fit(fd_fit.f(:), fd_fit.d(:), ftype, fopt);

% Output results.
varargout = {fitobject, gof, output};

end


