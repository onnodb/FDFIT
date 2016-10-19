function [d] = fOdijk(F, Lp, Lc, S, kT)
% Odijk Worm-like Chain
% NOTE: calculates d(F), not F(d)!
%
% INPUT:
% F = force (pN)
% Lp = persistence length (nm)
% Lc = contour length (um)
% S = stretching modulus (pN)
% kT = Boltzmann's constant times the temperature (optional; default:
%   4.11 pN nm)
%
% OUTPUT:
% d = extension (um)
%
% REFERENCES:
%   1. T. Odijk, Stiff Chains and Filaments under Tension, Macromolecules
%      28, 7016-7018 (1995).
%   2. M. D. Wang, H. Yin, R. Landick, J. Gelles, S. M. Block, Stretching
%      DNA with optical tweezers., Biophysical journal 72, 1335-46 (1997).

if nargin < 5
    kT = 4.11;
end

d = Lc .* (1 - 1./2*sqrt(kT./(F.*Lp)) + F./S);

% Disallow imaginary values
d(imag(d)~=0) = NaN;

end

