function [varargout] = fitfd(fd, varargin)
% FITFD Easy fitting of DNA force-extension curves.
%
% Based, by default, on least-squares fitting of the inverse of the eWLC/tWLC
% functions (i.e., expressing F as a function of d). Can also be called for
% fitting the regular way (with d as a function of F).
%
% SYNTAX:
% fitobject = fitfd(fd);
%   Most basic syntax; returns a cfit object (see MATLAB's 'fit' function).
% [fitobject, gof, output] = fitfd(fd);
%   Optionally, additionanl fit information can be returned (see MATLAB's
%   'fit' function).
% [fits] = fitfd(fdc);
%   Performs a fit for all FdData objects in an FdDataCollection.
%   A struct array with fit results is returned.
% fitfd(fd, 'key', value, ...);
%   Additional key-value pair arguments can be given; see below.
% fitfd(fd, 'flag', ...);
%   Additional flag arguments can be given; see below.
%
% INPUT:
% fd = FdData object. Note that when fitting the Odijk eWLC, the data is
%       automatically trimmed to only include forces up to 30 pN.
%       (See also the 'noTrim' flag argument below).
%
% KEY-VALUE PAIR ARGUMENTS:
% model = which model to use for fitting: one of the model classes
%       that descend from BasicFdFitModel (see directory Fitting/Models).
%
% FLAG ARGUMENTS:
% makePlot = automatically create a plot of the resulting fit.
% noTrim = when fitting the Odijk eWLC, do *not* automatically trim the data
%       to only include only forces up to 30 pN.
% skipFullFitInfo = only valid when fitting an FdDataCollection. When set, the
%       results struct array returned does not contain the fit information
%       fields "fd", "fitobject", "gof" and "output".
%       This is a bit of a hack, and seems to be necessary to prevent memory
%       corruption in older version of MATLAB when fitting an FdDataCollection.
%
% OUTPUT:
% fitobject = cfit object as returned by MATLAB's 'fit' function.
%       Depending on the model chosen (see the 'model' key-value pair argument),
%       the following fit coefficients/parameters are available:
%       - Lp: persistence length (nm)
%       - Lc: contour length (um)
%       - S: stretching modulus (pN)
%       - g0, g1: twist-stretch coupling constants (in pN*nm and nm,
%         respectively)
%       - F0: force offset (pN)
% gof = gof structure as returned by MATLAB's 'fit' function.
% output = output structure as returned by MATLAB's 'fit' function.
% fits = when fitting all items in an FdDataCollection, a struct array is
%       returned, where each item contains the fit parameters found.
%       Unless the 'skipFullFitInfo' flag is set, the following fields are also
%       available:
%       - fd: the FdData object fitted
%       - fitobject: see 'fitobject' output above
%       - gof: see 'fitobject' output above
%       - output: see 'fitobject' output above
%
% TODO:
% - Include FJC model
%
% SEE ALSO:
% fit, cfit, fOdijkInv, ftWLCInv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse input

defaultArgs = struct(...
                      'model',              [] ...
                    , 'makePlot',           false ...
                    , 'noTrim',             false ...
                    , 'skipFullFitInfo',    false ...
                    );

args = parseArgs(varargin, defaultArgs, {'makePlot','noTrim','skipFullFitInfo'});

if isempty(args.model) || ~isa(args.model, 'BasicFdFitModel')
    error('fitfd:InvalidFitModel', 'Invalid fit model.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process FdDataCollection if given

if isa(fd, 'FdDataCollection')
    % Pre-allocate results.
    res = struct('fd', []);
    res(fd.length).fd = [];

    % Do the fits.
    parfor i = 1:fd.length
        [f,g,o] = fitfd(fd.items{i}, varargin{:}); %#ok

        % Store fit parameters directly.
        cn = coeffnames(f);
        for k = 1:length(cn)
            res(i).(cn{k}) = f.(cn{k});
        end

        % Store more detailed fit information, too, if requested.
        if ~args.skipFullFitInfo %#ok
            res(i).fd        = fd.items{i};
            res(i).fitobject = f;
            res(i).gof       = g;
            res(i).output    = o;
        end
    end
    if args.skipFullFitInfo
        res = rmfield(res, 'fd');
    end
    varargout = {res};
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do fit

% Check if data is OK for fitting (and possibly filter data).
fd_fit = args.model.validateData(fd);

% Do actual fit.
ftype = args.model.getFitType();
fopt  = args.model.getFitOptions();
fopt.MaxIter    = 120;

switch args.model.dependentVariable
    case 'F'
        [fitobject, gof, output] = fit(fd_fit.d(:), fd_fit.f(:), ftype, fopt);
    case 'd'
        [fitobject, gof, output] = fit(fd_fit.f(:), fd_fit.d(:), ftype, fopt);
    otherwise
        error('Invalid dependent variable for model object.');
end

% Make plot, if requested.
if args.makePlot
    plotfdfit(fd, fitobject);
end

% Output results.
varargout = {fitobject, gof, output};

end

