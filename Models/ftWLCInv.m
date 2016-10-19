function [F] = ftWLCInv(d, Lp, Lc, S, C, g0, g1, Fc, kT)
% Twistable Worm-like Chain (inverse)
%
% INPUT:
% d = extension (can be a single, scalar value, or a vector) (um)
% Lp = persistence length (nm)
% Lc = contour length (um)
% S = stretching modulus (pN)
% C = twist rigidity (pN nm^2)
% g0, g1 = twist-stretch coupling (g0: pN Nm; g1: nm)
% Fc = critical force for twist-stretch coupling (pN)
% kT = Boltzmann's constant times the temperature (optional; default:
%   4.11 pN nm)
%
% OUTPUT:
% F = force (pN)
%
% NOTES:
% For a single 'd' value, the 'ftWLC' function is inverted numerically, using
% the built-in 'fminbnd' function. This is relatively slow, so if you're, e.g.,
% plotting a graph, you should call this function with a vector of 'd'
% values. In that case, the inversion is performed using a spline-interpolated
% lookup table.
%
% SEE ALSO:
% ftWLC

if nargin < 9
    kT = 4.11;
end

Fmax = (-g0+sqrt(S*C))/g1;      % above Fmax, the tWLC function becomes rubbish

if isscalar(d)
    % For a single distance value, we just invert the tWLC numerically using
    % 'fminbnd'
    F = fminbnd(@optimFun, 0, Fmax);
else
    % For a vector, we use a more optimized version, using an interpolated
    % lookup table
    dF = 1e-2;                  % spacing for lookup table [1]
    f_lookup = dF:dF:Fmax;
    if isempty(f_lookup)
        % Invalid parameters: impossible to generate model curve for
        % interpolation
        F = zeros(size(d));
    else
        d_lookup = ftWLC(f_lookup, Lp, Lc, S, C, g0, g1, Fc, kT);
        F = interp1(d_lookup, f_lookup, d, 'spline');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [y] = optimFun(F_trial)
        y = ftWLC(F_trial, Lp, Lc, S, C, g0, g1, Fc, kT) - d;
        y(isnan(y) | (imag(y) ~= 0)) = +Inf;
        y = abs(y);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOOTNOTES:
% [1] a 'dF' value of 1e-2 is sufficient to have only an error of order
%     1e-5 compared to the 'fminbnd' strategy, for typical parameter
%     values.

end
