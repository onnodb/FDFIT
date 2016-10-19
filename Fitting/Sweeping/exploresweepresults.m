function exploresweepresults(fd, sweepResults)
% EXPLORESWEEPRESULTS Convenience function for exploring the results of a 1D sweep.
%
% SYNTAX:
% exploresweepresults(fd, sweepResults);
%
% INPUT:
% fd = an FdData object.
% sweepResults = result of a call to the "sweepfdfits" function, on the FdData
%       object in "fd".
%
% EXAMPLE:
% >> sr = sweepfdfits(fd, ...);
% >> exploresweepresults(fd, sr);
%
% SEE ALSO:
% sweepfdfits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

switch sweepResults.options.sweepType
    case {'left','right'}
        % ok
    otherwise
        error('plotsweepresults can only be used on the results of a one-boundary sweep.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use 'explorecell' to do the actual work

explorecell(...
    sweepResults.params.all, ...
    @(ax,item,idx) plotfdfit(...
                             ax, fd, item, ...
                             'model', sweepResults.options.model, ...
                             'Lc', sweepResults.options.Lc, ...
                             'highlightSubset', fd.subset(sweepResults.options.sweepBoundaryType, ...
                                                          sweepResults.sweeps(idx).dataRange) ...
                             ) ...
    );

end

