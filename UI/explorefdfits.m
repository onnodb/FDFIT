function explorefdfits(fits, varargin)
% EXPLOREFDFITS Convenience function for inspecting multiple 'fitfd' results.
%
% EXAMPLE:
% fits = fitfd(fdc, 'model', 'twlc-f0');
% explorefdfits(fits, 'model', 'twlc-f0');
%   ...where 'fdc' is an FdDataCollection object.
%
% NOTE:
% If the 'fits' argument was obtained by a call to 'fitfd' with the
% 'skipFullFitInfo' flag set, this function will fail.
%
% SEE ALSO:
% fitfd, explorecell

explorecell(...
    length(fits), ...
    @(ax,item,idx) plotfdfit(ax, fits(idx).fd, fits(idx).fitobject, varargin{:}) ...
    );

end

