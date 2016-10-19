function [d] = ftWLC(F, Lp, Lc, S, C, g0, g1, Fc, kT)
% Twistable Worm-like Chain
% NOTE: calculates d(F), not F(d)!
%
% INPUT:
% F = force (pN)
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
% d = extension (um)
%
% REFERENCES:
%   1. P. Gross et al., Quantifying how DNA stretches, melts and changes
%      twist under tension, Nature Physics 7, 731-736 (2011).

if nargin < 9
    kT = 4.11;
end

if isscalar(F)
    if F < Fc
        g = g0 + g1 * Fc;
    else
        g = g0 + g1 * F;
    end
else
    g = zeros(size(F));
    g(F < Fc)  = g0 + g1 * Fc;
    g(F >= Fc) = g0 + g1 .* F(F >= Fc);
end

d = Lc .* (1 - 1./2*sqrt(kT./(F.*Lp)) + (C ./ (-g.^2+S.*C)) .* F );

end

