function [d] = fFJC(F, Lp, Lc, S, kT)
% Freely-Jointed Chain
% NOTE: calculates d(F), not F(d)!
%
% INPUT:
% F = force (pN)
% Lp = persistence length (nm)
% Lc = contour length (um)
% S = elastic modulus (pN)
% kT = Boltzmann's constant times the temperature (optional; default:
%   4.11 pN nm)
%
% OUTPUT:
% d = extension (um)
%
% REFERENCES:
%   1. S. B. Smith, Y. Cui, C. Bustamante, Overstretching B-DNA: The
%      Elastic Response of Individual Double-Stranded and Single-Stranded
%      DNA Molecules, Science 271, 795-799 (1996).
%   2. M. D. Wang, H. Yin, R. Landick, J. Gelles, S. M. Block, Stretching
%      DNA with optical tweezers., Biophysical journal 72, 1335-46 (1997).

if nargin < 5
    kT = 4.11;
end

d = Lc ...
    .* ...
        ( ...
            coth(2.*F.*Lp./kT)  ...
            - kT./(2.*F.*Lp)    ...
        ) ...
    .* ...
    (1 + F./S);

end
