function [d] = fOdijk_f0(F, Lp, Lc, S, F0, kT)

if nargin < 6
    kT = 4.11;
end

d = fOdijk(F-F0, Lp, Lc, S, kT);

end
