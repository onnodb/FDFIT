function [F] = fOdijkInv_d0(d, Lp, Lc, S, d0, kT)
% Odijk Worm-like Chain (inverse) with distance offset
%
% INPUT:
% d = extension (um)
% Lp = persistence length (nm)
% Lc = contour length (um)
% S = stretching modulus (pN)
% d0 = distance offset (um)
% kT = Boltzmann's constant times the temperature (optional; default:
%   4.11 pN nm)
%
% OUTPUT:
% F = force (pN)
%
% SEE ALSO:
% fOdijkInv

if nargin < 6
    kT = 4.11;
end

F = fOdijkInv(d+d0, Lp, Lc, S, kT);

end
